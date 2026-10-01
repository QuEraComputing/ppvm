// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`CudaPauliSum`]: a Pauli sum stored on the GPU and propagated with
//! cuPauliProp.
//!
//! cuPauliProp applies operators out of place, so the state is two
//! pre-allocated device expansions used alternately as input and output.
//! Every apply deduplicates and truncates (`|c| <= cutoff`), matching ppvm's
//! truncate-after-every-gate.

use std::collections::HashMap;
use std::ffi::{CStr, c_void};
use std::ptr;

use ppvm_cupauliprop_sys::*;

use crate::Propagator;
use crate::encode::{num_words, pauli_channel_probs, push_packed_term};

fn check(status: cupaulipropStatus_t, call: &str) {
    if status != CUPAULIPROP_STATUS_SUCCESS {
        let msg = unsafe { CStr::from_ptr(cupaulipropGetErrorString(status)) };
        panic!(
            "{call} failed with status {status}: {}",
            msg.to_string_lossy()
        );
    }
}

fn check_cuda(err: cudaError_t, call: &str) {
    if err != 0 {
        let msg = unsafe { CStr::from_ptr(cudaGetErrorString(err)) };
        panic!(
            "{call} failed with CUDA error {err}: {}",
            msg.to_string_lossy()
        );
    }
}

fn device_alloc(bytes: usize) -> *mut c_void {
    let mut p = ptr::null_mut();
    check_cuda(unsafe { cudaMalloc(&mut p, bytes.max(16)) }, "cudaMalloc");
    p
}

/// `cupaulipropGetVersion` as `major.minor.patch`.
pub fn library_version() -> String {
    let v = unsafe { cupaulipropGetVersion() };
    format!("{}.{}.{}", v / 10000, (v / 100) % 100, v % 100)
}

/// `(free, total)` device memory in bytes.
pub fn device_memory() -> (usize, usize) {
    let (mut free, mut total) = (0, 0);
    check_cuda(
        unsafe { cudaMemGetInfo(&mut free, &mut total) },
        "cudaMemGetInfo",
    );
    (free, total)
}

/// Device memory sizing. Unset fields are derived from free device memory.
#[derive(Debug, Clone, Copy, Default)]
pub struct CudaOptions {
    /// Term capacity of each of the two expansions.
    pub capacity: Option<usize>,
    /// Workspace buffer size in bytes.
    pub workspace_bytes: Option<usize>,
}

impl CudaOptions {
    /// Resolve to `(capacity, workspace_bytes)` for `n_qubits`. By default
    /// 90% of free memory is split evenly between the workspace and the two
    /// expansions.
    fn resolve(self, n_qubits: usize) -> (usize, usize) {
        let term_bytes = 16 * num_words(n_qubits) + 8;
        let budget = device_memory().0 / 10 * 9;
        let workspace = self.workspace_bytes.unwrap_or(budget / 3);
        let capacity = self
            .capacity
            .unwrap_or_else(|| budget.saturating_sub(workspace) / 2 / term_bytes);
        (capacity, workspace)
    }
}

/// Cache key for created operators: kind, qubits and parameter bit patterns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum OpKey {
    Rotation {
        qubits: [i32; 2],
        paulis: [cupaulipropPauliKind_t; 2],
        n: usize,
        angle: u64,
    },
    PauliChannel {
        qubit: i32,
        probs: [u64; 3],
    },
}

struct DeviceExpansion {
    expansion: cupaulipropPauliExpansion_t,
    xz: *mut c_void,
    coef: *mut c_void,
}

/// A Pauli sum on the GPU, propagated backwards (Heisenberg picture) like
/// ppvm's `PauliSum`. Gate methods take the same arguments as ppvm's.
pub struct CudaPauliSum {
    handle: cupaulipropHandle_t,
    workspace: cupaulipropWorkspaceDescriptor_t,
    workspace_buf: *mut c_void,
    workspace_bytes: usize,
    expansions: [DeviceExpansion; 2],
    /// Index of the expansion holding the current state.
    current: usize,
    n_qubits: usize,
    capacity: usize,
    cutoff: f64,
    operators: HashMap<OpKey, cupaulipropQuantumOperator_t>,
}

impl CudaPauliSum {
    /// Upload dense Pauli strings (all coefficients 1) as the initial
    /// observable, truncating at `|c| <= cutoff` on every subsequent gate.
    pub fn new(n_qubits: usize, terms: &[String], cutoff: f64, options: CudaOptions) -> Self {
        let words = num_words(n_qubits);
        let mut lib_words = 0;
        check(
            unsafe { cupaulipropGetNumPackedIntegers(n_qubits as i32, &mut lib_words) },
            "cupaulipropGetNumPackedIntegers",
        );
        assert_eq!(lib_words as usize, words, "packed-word count mismatch");

        let (capacity, workspace_bytes) = options.resolve(n_qubits);
        assert!(
            terms.len() <= capacity,
            "{} initial terms exceed capacity {capacity}",
            terms.len()
        );

        let mut handle = ptr::null_mut();
        check(
            unsafe { cupaulipropCreate(&mut handle) },
            "cupaulipropCreate",
        );
        let mut workspace = ptr::null_mut();
        check(
            unsafe { cupaulipropCreateWorkspaceDescriptor(handle, &mut workspace) },
            "cupaulipropCreateWorkspaceDescriptor",
        );
        let workspace_buf = device_alloc(workspace_bytes);

        let mut packed = Vec::with_capacity(terms.len() * 2 * words);
        for term in terms {
            assert_eq!(term.len(), n_qubits, "term {term:?} has wrong length");
            push_packed_term(&mut packed, term, words);
        }
        let coefs = vec![1.0f64; terms.len()];

        let xz_bytes = capacity * 2 * words * 8;
        let coef_bytes = capacity * 8;
        let make = |n_terms: usize| {
            let xz = device_alloc(xz_bytes);
            let coef = device_alloc(coef_bytes);
            let mut expansion = ptr::null_mut();
            check(
                unsafe {
                    cupaulipropCreatePauliExpansion(
                        handle,
                        n_qubits as i32,
                        xz,
                        xz_bytes as i64,
                        coef,
                        coef_bytes as i64,
                        CUDA_R_64F,
                        n_terms as i64,
                        CUPAULIPROP_SORT_ORDER_NONE,
                        0,
                        &mut expansion,
                    )
                },
                "cupaulipropCreatePauliExpansion",
            );
            DeviceExpansion {
                expansion,
                xz,
                coef,
            }
        };

        let first = make(terms.len());
        unsafe {
            check_cuda(
                cudaMemcpy(
                    first.xz,
                    packed.as_ptr().cast(),
                    packed.len() * 8,
                    CUDA_MEMCPY_HOST_TO_DEVICE,
                ),
                "cudaMemcpy",
            );
            check_cuda(
                cudaMemcpy(
                    first.coef,
                    coefs.as_ptr().cast(),
                    coefs.len() * 8,
                    CUDA_MEMCPY_HOST_TO_DEVICE,
                ),
                "cudaMemcpy",
            );
        }
        let second = make(0);

        Self {
            handle,
            workspace,
            workspace_buf,
            workspace_bytes,
            expansions: [first, second],
            current: 0,
            n_qubits,
            capacity,
            cutoff,
            operators: HashMap::new(),
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Term capacity of each expansion.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Wait for all queued GPU work.
    pub fn synchronize(&self) {
        check_cuda(unsafe { cudaDeviceSynchronize() }, "cudaDeviceSynchronize");
    }

    fn num_terms(&self) -> usize {
        let mut n = 0;
        check(
            unsafe {
                cupaulipropPauliExpansionGetNumTerms(
                    self.handle,
                    self.expansions[self.current].expansion,
                    &mut n,
                )
            },
            "cupaulipropPauliExpansionGetNumTerms",
        );
        n as usize
    }

    /// View over all terms of the current expansion; destroy after use.
    fn current_view(&self) -> cupaulipropPauliExpansionView_t {
        let mut view = ptr::null_mut();
        check(
            unsafe {
                cupaulipropPauliExpansionGetContiguousRange(
                    self.handle,
                    self.expansions[self.current].expansion,
                    0,
                    self.num_terms() as i64,
                    &mut view,
                )
            },
            "cupaulipropPauliExpansionGetContiguousRange",
        );
        view
    }

    /// Check the size a `Prepare*` call recorded and (re)attach the buffer,
    /// since `Prepare*` detaches it.
    fn attach_workspace(&self, what: &str) {
        let mut required = 0;
        check(
            unsafe {
                cupaulipropWorkspaceGetMemorySize(
                    self.handle,
                    self.workspace,
                    CUPAULIPROP_MEMSPACE_DEVICE,
                    CUPAULIPROP_WORKSPACE_SCRATCH,
                    &mut required,
                )
            },
            "cupaulipropWorkspaceGetMemorySize",
        );
        assert!(
            required as usize <= self.workspace_bytes,
            "{what} needs {required} workspace bytes, have {}",
            self.workspace_bytes
        );
        check(
            unsafe {
                cupaulipropWorkspaceSetMemory(
                    self.handle,
                    self.workspace,
                    CUPAULIPROP_MEMSPACE_DEVICE,
                    CUPAULIPROP_WORKSPACE_SCRATCH,
                    self.workspace_buf,
                    self.workspace_bytes as i64,
                )
            },
            "cupaulipropWorkspaceSetMemory",
        );
    }

    fn operator(&mut self, key: OpKey) -> cupaulipropQuantumOperator_t {
        let handle = self.handle;
        *self.operators.entry(key).or_insert_with(|| {
            let mut op = ptr::null_mut();
            let status = match key {
                OpKey::Rotation {
                    qubits,
                    paulis,
                    n,
                    angle,
                } => unsafe {
                    cupaulipropCreatePauliRotationGateOperator(
                        handle,
                        f64::from_bits(angle),
                        n as i32,
                        qubits.as_ptr(),
                        paulis.as_ptr(),
                        &mut op,
                    )
                },
                OpKey::PauliChannel { qubit, probs } => {
                    let probs = pauli_channel_probs(probs.map(f64::from_bits));
                    unsafe {
                        cupaulipropCreatePauliNoiseChannelOperator(
                            handle,
                            1,
                            &qubit,
                            probs.as_ptr(),
                            &mut op,
                        )
                    }
                }
            };
            check(status, "cupaulipropCreate*Operator");
            op
        })
    }

    /// Apply the adjoint of `key` (Heisenberg picture), then deduplicate and
    /// truncate into the other expansion.
    fn apply(&mut self, key: OpKey) {
        let op = self.operator(key);
        let mut params = cupaulipropCoefficientTruncationParams_t {
            cutoff: self.cutoff,
        };
        let strategy = cupaulipropTruncationStrategy_t {
            strategy: CUPAULIPROP_TRUNCATION_STRATEGY_COEFFICIENT_BASED,
            paramStruct: (&raw mut params).cast(),
        };
        let view = self.current_view();

        let (mut xz_needed, mut coef_needed) = (0, 0);
        check(
            unsafe {
                cupaulipropPauliExpansionViewPrepareOperatorApplication(
                    self.handle,
                    view,
                    op,
                    CUPAULIPROP_SORT_ORDER_INTERNAL,
                    0,
                    1,
                    &strategy,
                    self.workspace_bytes as i64,
                    &mut xz_needed,
                    &mut coef_needed,
                    self.workspace,
                )
            },
            "cupaulipropPauliExpansionViewPrepareOperatorApplication",
        );
        let needed_terms = (coef_needed as usize).div_ceil(8);
        assert!(
            needed_terms <= self.capacity,
            "operator application needs capacity {needed_terms} terms, have {}",
            self.capacity
        );
        self.attach_workspace("operator application");

        let out = self.expansions[1 - self.current].expansion;
        check(
            unsafe {
                cupaulipropPauliExpansionViewComputeOperatorApplication(
                    self.handle,
                    view,
                    out,
                    op,
                    1,
                    CUPAULIPROP_SORT_ORDER_INTERNAL,
                    0,
                    1,
                    &strategy,
                    self.workspace,
                    ptr::null_mut(),
                )
            },
            "cupaulipropPauliExpansionViewComputeOperatorApplication",
        );
        check(
            unsafe { cupaulipropDestroyPauliExpansionView(view) },
            "cupaulipropDestroyPauliExpansionView",
        );
        self.current = 1 - self.current;
    }

    fn rotation(&mut self, qubits: &[usize], paulis: &[cupaulipropPauliKind_t], theta: f64) {
        let mut key_qubits = [0; 2];
        let mut key_paulis = [CUPAULIPROP_PAULI_I; 2];
        for (i, (&q, &p)) in qubits.iter().zip(paulis).enumerate() {
            key_qubits[i] = q as i32;
            key_paulis[i] = p;
        }
        self.apply(OpKey::Rotation {
            qubits: key_qubits,
            paulis: key_paulis,
            n: qubits.len(),
            angle: theta.to_bits(),
        });
    }

    pub fn rx(&mut self, q: usize, theta: f64) {
        self.rotation(&[q], &[CUPAULIPROP_PAULI_X], theta);
    }

    pub fn ry(&mut self, q: usize, theta: f64) {
        self.rotation(&[q], &[CUPAULIPROP_PAULI_Y], theta);
    }

    pub fn rz(&mut self, q: usize, theta: f64) {
        self.rotation(&[q], &[CUPAULIPROP_PAULI_Z], theta);
    }

    pub fn rxx(&mut self, a: usize, b: usize, theta: f64) {
        self.rotation(&[a, b], &[CUPAULIPROP_PAULI_X; 2], theta);
    }

    pub fn ryy(&mut self, a: usize, b: usize, theta: f64) {
        self.rotation(&[a, b], &[CUPAULIPROP_PAULI_Y; 2], theta);
    }

    pub fn rzz(&mut self, a: usize, b: usize, theta: f64) {
        self.rotation(&[a, b], &[CUPAULIPROP_PAULI_Z; 2], theta);
    }

    /// Single-qubit Pauli channel with probabilities `[p_x, p_y, p_z]`.
    pub fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        self.apply(OpKey::PauliChannel {
            qubit: q as i32,
            probs: p.map(f64::to_bits),
        });
    }

    /// `Tr(O |0…0⟩⟨0…0|)`, the sum of coefficients of `{I, Z}`-only strings.
    pub fn overlap_with_zero(&self) -> f64 {
        let view = self.current_view();
        check(
            unsafe {
                cupaulipropPauliExpansionViewPrepareTraceWithZeroState(
                    self.handle,
                    view,
                    self.workspace_bytes as i64,
                    self.workspace,
                )
            },
            "cupaulipropPauliExpansionViewPrepareTraceWithZeroState",
        );
        self.attach_workspace("trace");
        let (mut significand, mut exponent) = (0.0f64, 0.0f64);
        check(
            unsafe {
                cupaulipropPauliExpansionViewComputeTraceWithZeroState(
                    self.handle,
                    view,
                    (&raw mut significand).cast(),
                    &mut exponent,
                    self.workspace,
                    ptr::null_mut(),
                )
            },
            "cupaulipropPauliExpansionViewComputeTraceWithZeroState",
        );
        self.synchronize();
        check(
            unsafe { cupaulipropDestroyPauliExpansionView(view) },
            "cupaulipropDestroyPauliExpansionView",
        );
        significand * exponent.exp2()
    }

    pub fn len(&self) -> usize {
        self.num_terms()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Drop for CudaPauliSum {
    fn drop(&mut self) {
        // Ignore statuses: panicking in drop would abort.
        unsafe {
            for op in self.operators.values() {
                cupaulipropDestroyOperator(*op);
            }
            for e in &self.expansions {
                cupaulipropDestroyPauliExpansion(e.expansion);
                cudaFree(e.xz);
                cudaFree(e.coef);
            }
            cupaulipropDestroyWorkspaceDescriptor(self.workspace);
            cudaFree(self.workspace_buf);
            cupaulipropDestroy(self.handle);
        }
    }
}

impl Propagator for CudaPauliSum {
    fn rx(&mut self, q: usize, theta: f64) {
        CudaPauliSum::rx(self, q, theta);
    }

    fn rxx(&mut self, a: usize, b: usize, theta: f64) {
        CudaPauliSum::rxx(self, a, b, theta);
    }

    fn ryy(&mut self, a: usize, b: usize, theta: f64) {
        CudaPauliSum::ryy(self, a, b, theta);
    }

    fn rzz(&mut self, a: usize, b: usize, theta: f64) {
        CudaPauliSum::rzz(self, a, b, theta);
    }

    fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        CudaPauliSum::pauli_error(self, q, p);
    }

    fn overlap_with_zero(&self) -> f64 {
        CudaPauliSum::overlap_with_zero(self)
    }

    fn len(&self) -> usize {
        CudaPauliSum::len(self)
    }
}
