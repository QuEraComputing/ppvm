// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Raw FFI declarations for the subset of NVIDIA cuPauliProp (cuQuantum,
//! `cupauliprop.h` version 0.5) and the CUDA runtime used by
//! `ppvm-cupauliprop`.
//!
//! Names follow the C header so declarations can be checked against it line
//! by line. C enums are declared as `i32` aliases plus constants. Nothing is
//! linked unless the `cuda` feature is enabled (see `build.rs`).
#![allow(non_camel_case_types, non_snake_case)]

use std::ffi::{c_char, c_void};

// ---- cuPauliProp opaque handles ----
pub type cupaulipropHandle_t = *mut c_void;
pub type cupaulipropWorkspaceDescriptor_t = *mut c_void;
pub type cupaulipropPauliExpansion_t = *mut c_void;
pub type cupaulipropPauliExpansionView_t = *mut c_void;
pub type cupaulipropQuantumOperator_t = *mut c_void;
pub type cupaulipropPackedIntegerType_t = u64;

// ---- CUDA runtime types ----
pub type cudaStream_t = *mut c_void;
pub type cudaError_t = i32;
pub type cudaDataType_t = i32;
pub const CUDA_R_64F: cudaDataType_t = 1;

pub type cudaMemcpyKind = i32;
pub const CUDA_MEMCPY_HOST_TO_DEVICE: cudaMemcpyKind = 1;
pub const CUDA_MEMCPY_DEVICE_TO_HOST: cudaMemcpyKind = 2;

// ---- cuPauliProp enums ----
pub type cupaulipropStatus_t = i32;
pub const CUPAULIPROP_STATUS_SUCCESS: cupaulipropStatus_t = 0;

pub type cupaulipropMemspace_t = i32;
pub const CUPAULIPROP_MEMSPACE_DEVICE: cupaulipropMemspace_t = 0;

pub type cupaulipropWorkspaceKind_t = i32;
pub const CUPAULIPROP_WORKSPACE_SCRATCH: cupaulipropWorkspaceKind_t = 0;

pub type cupaulipropTruncationStrategyKind_t = i32;
pub const CUPAULIPROP_TRUNCATION_STRATEGY_COEFFICIENT_BASED: cupaulipropTruncationStrategyKind_t =
    0;
pub const CUPAULIPROP_TRUNCATION_STRATEGY_PAULI_WEIGHT_BASED: cupaulipropTruncationStrategyKind_t =
    1;

pub type cupaulipropSortOrder_t = i32;
pub const CUPAULIPROP_SORT_ORDER_NONE: cupaulipropSortOrder_t = 0;
pub const CUPAULIPROP_SORT_ORDER_INTERNAL: cupaulipropSortOrder_t = 1;

pub type cupaulipropPauliKind_t = i32;
pub const CUPAULIPROP_PAULI_I: cupaulipropPauliKind_t = 0;
pub const CUPAULIPROP_PAULI_X: cupaulipropPauliKind_t = 1;
pub const CUPAULIPROP_PAULI_Y: cupaulipropPauliKind_t = 2;
pub const CUPAULIPROP_PAULI_Z: cupaulipropPauliKind_t = 3;

// ---- cuPauliProp structs ----
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct cupaulipropTruncationStrategy_t {
    pub strategy: cupaulipropTruncationStrategyKind_t,
    pub paramStruct: *mut c_void,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct cupaulipropCoefficientTruncationParams_t {
    pub cutoff: f64,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct cupaulipropPauliWeightTruncationParams_t {
    pub cutoff: i32,
}

unsafe extern "C" {
    // ---- library management ----
    pub fn cupaulipropGetVersion() -> usize;
    pub fn cupaulipropGetErrorString(error: cupaulipropStatus_t) -> *const c_char;
    pub fn cupaulipropGetNumPackedIntegers(
        numQubits: i32,
        numPackedIntegers: *mut i32,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropCreate(handle: *mut cupaulipropHandle_t) -> cupaulipropStatus_t;
    pub fn cupaulipropDestroy(handle: cupaulipropHandle_t) -> cupaulipropStatus_t;

    // ---- workspace ----
    pub fn cupaulipropCreateWorkspaceDescriptor(
        handle: cupaulipropHandle_t,
        workspaceDesc: *mut cupaulipropWorkspaceDescriptor_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropDestroyWorkspaceDescriptor(
        workspaceDesc: cupaulipropWorkspaceDescriptor_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropWorkspaceGetMemorySize(
        handle: cupaulipropHandle_t,
        workspaceDesc: cupaulipropWorkspaceDescriptor_t,
        memSpace: cupaulipropMemspace_t,
        workspaceKind: cupaulipropWorkspaceKind_t,
        memoryBufferSize: *mut i64,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropWorkspaceSetMemory(
        handle: cupaulipropHandle_t,
        workspaceDesc: cupaulipropWorkspaceDescriptor_t,
        memSpace: cupaulipropMemspace_t,
        workspaceKind: cupaulipropWorkspaceKind_t,
        memoryBuffer: *mut c_void,
        memoryBufferSize: i64,
    ) -> cupaulipropStatus_t;

    // ---- Pauli expansions and views ----
    pub fn cupaulipropCreatePauliExpansion(
        handle: cupaulipropHandle_t,
        numQubits: i32,
        xzBitsBuffer: *mut c_void,
        xzBitsBufferSize: i64,
        coefBuffer: *mut c_void,
        coefBufferSize: i64,
        dataType: cudaDataType_t,
        numLocalTerms: i64,
        sortOrder: cupaulipropSortOrder_t,
        hasDuplicates: i32,
        pauliExpansion: *mut cupaulipropPauliExpansion_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropDestroyPauliExpansion(
        pauliExpansion: cupaulipropPauliExpansion_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropPauliExpansionGetNumTerms(
        handle: cupaulipropHandle_t,
        pauliExpansion: cupaulipropPauliExpansion_t,
        numLocalTerms: *mut i64,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropPauliExpansionGetContiguousRange(
        handle: cupaulipropHandle_t,
        pauliExpansion: cupaulipropPauliExpansion_t,
        startIndex: i64,
        endIndex: i64,
        view: *mut cupaulipropPauliExpansionView_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropDestroyPauliExpansionView(
        view: cupaulipropPauliExpansionView_t,
    ) -> cupaulipropStatus_t;

    // ---- trace with |0…0⟩ ----
    pub fn cupaulipropPauliExpansionViewPrepareTraceWithZeroState(
        handle: cupaulipropHandle_t,
        view: cupaulipropPauliExpansionView_t,
        maxWorkspaceDeviceSize: i64,
        workspace: cupaulipropWorkspaceDescriptor_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropPauliExpansionViewComputeTraceWithZeroState(
        handle: cupaulipropHandle_t,
        view: cupaulipropPauliExpansionView_t,
        traceSignificand: *mut c_void,
        traceExponent: *mut f64,
        workspace: cupaulipropWorkspaceDescriptor_t,
        stream: cudaStream_t,
    ) -> cupaulipropStatus_t;

    // ---- operator application ----
    pub fn cupaulipropPauliExpansionViewPrepareOperatorApplication(
        handle: cupaulipropHandle_t,
        viewIn: cupaulipropPauliExpansionView_t,
        quantumOperator: cupaulipropQuantumOperator_t,
        sortOrder: cupaulipropSortOrder_t,
        keepDuplicates: i32,
        numTruncationStrategies: i32,
        truncationStrategies: *const cupaulipropTruncationStrategy_t,
        maxWorkspaceDeviceSize: i64,
        requiredXZBitsBufferSize: *mut i64,
        requiredCoefBufferSize: *mut i64,
        workspace: cupaulipropWorkspaceDescriptor_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropPauliExpansionViewComputeOperatorApplication(
        handle: cupaulipropHandle_t,
        viewIn: cupaulipropPauliExpansionView_t,
        expansionOut: cupaulipropPauliExpansion_t,
        quantumOperator: cupaulipropQuantumOperator_t,
        adjoint: i32,
        sortOrder: cupaulipropSortOrder_t,
        keepDuplicates: i32,
        numTruncationStrategies: i32,
        truncationStrategies: *const cupaulipropTruncationStrategy_t,
        workspace: cupaulipropWorkspaceDescriptor_t,
        stream: cudaStream_t,
    ) -> cupaulipropStatus_t;

    // ---- quantum operators ----
    pub fn cupaulipropCreatePauliRotationGateOperator(
        handle: cupaulipropHandle_t,
        angle: f64,
        numQubits: i32,
        qubitIndices: *const i32,
        paulis: *const cupaulipropPauliKind_t,
        oper: *mut cupaulipropQuantumOperator_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropCreatePauliNoiseChannelOperator(
        handle: cupaulipropHandle_t,
        numQubits: i32,
        qubitIndices: *const i32,
        probabilities: *const f64,
        oper: *mut cupaulipropQuantumOperator_t,
    ) -> cupaulipropStatus_t;
    pub fn cupaulipropDestroyOperator(oper: cupaulipropQuantumOperator_t) -> cupaulipropStatus_t;

    // ---- CUDA runtime ----
    pub fn cudaMalloc(devPtr: *mut *mut c_void, size: usize) -> cudaError_t;
    pub fn cudaFree(devPtr: *mut c_void) -> cudaError_t;
    pub fn cudaMemcpy(
        dst: *mut c_void,
        src: *const c_void,
        count: usize,
        kind: cudaMemcpyKind,
    ) -> cudaError_t;
    pub fn cudaMemGetInfo(free: *mut usize, total: *mut usize) -> cudaError_t;
    pub fn cudaDeviceSynchronize() -> cudaError_t;
    pub fn cudaGetErrorString(error: cudaError_t) -> *const c_char;
}
