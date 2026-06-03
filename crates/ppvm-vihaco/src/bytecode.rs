//! Binary bytecode (`.ssb`) round-trip for PPVM modules.
//!
//! Serializes a resolved [`Module`] to a little-endian container — header
//! (magic + version + device info) → strings section → code section — and
//! reads it back. The per-instruction codec is the existing `WriteBytes` /
//! `FromBytes` generated for [`PPVMInstruction`]; this module only frames the
//! container around it.

use std::io::{Read, Write};

use vihaco::instruction::{FromBytes, WriteBytes};

use crate::PPVMModule;
use crate::composite::{PPVM_MAGIC, PPVMDeviceInfo, PPVMInstruction};

/// Current `.ssb` format version. The reader rejects any other version.
pub const PPVM_BYTECODE_VERSION: u16 = 1;

/// Byte length of the fixed v1 header (magic 4, version 2, header_size 4,
/// n_qubits 4, coefficient_threshold 8) and the offset where the strings
/// section begins.
const HEADER_SIZE: u32 = 4 + 2 + 4 + 4 + 8;

/// Serialize a resolved module to the v1 `.ssb` byte stream.
pub fn write_module<W: Write>(module: &PPVMModule, w: &mut W) -> eyre::Result<()> {
    // v1 serializes only code, strings, and device info. Refuse to silently
    // drop any table a future feature might populate.
    let populated = if !module.functions.is_empty() {
        Some("functions")
    } else if !module.labels.is_empty() {
        Some("labels")
    } else if !module.constants.is_empty() {
        Some("constants")
    } else if !module.source_symbols.is_empty() {
        Some("source_symbols")
    } else if module.main_function.is_some() {
        Some("main_function")
    } else if module.file != 0 {
        Some("file")
    } else {
        None
    };
    if let Some(table) = populated {
        return Err(eyre::eyre!(
            "bytecode v1 cannot represent a populated `{table}`"
        ));
    }

    let info = &module.extra;
    let n_qubits = u32::try_from(info.n_qubits)
        .map_err(|_| eyre::eyre!("n_qubits {} does not fit in u32", info.n_qubits))?;

    // Header.
    w.write_all(&PPVM_MAGIC.to_le_bytes())?;
    w.write_all(&PPVM_BYTECODE_VERSION.to_le_bytes())?;
    w.write_all(&HEADER_SIZE.to_le_bytes())?;
    w.write_all(&n_qubits.to_le_bytes())?;
    w.write_all(&info.coefficient_threshold.to_le_bytes())?;

    // Strings section: count, then each entry as len-prefixed UTF-8.
    let string_count =
        u32::try_from(module.strings.len()).map_err(|_| eyre::eyre!("string count exceeds u32"))?;
    w.write_all(&string_count.to_le_bytes())?;
    for s in &module.strings {
        let len = u32::try_from(s.len()).map_err(|_| eyre::eyre!("string length exceeds u32"))?;
        w.write_all(&len.to_le_bytes())?;
        w.write_all(s.as_bytes())?;
    }

    // Code section: count, then each instruction's fixed-width frame.
    let code_count =
        u32::try_from(module.code.len()).map_err(|_| eyre::eyre!("code length exceeds u32"))?;
    w.write_all(&code_count.to_le_bytes())?;
    for inst in &module.code {
        inst.write_bytes(w)?;
    }

    Ok(())
}

/// Reconstruct a module from a v1 `.ssb` byte stream.
pub fn read_module<R: Read>(r: &mut R) -> eyre::Result<PPVMModule> {
    // Header.
    let magic = read_u32(r)?;
    if magic != PPVM_MAGIC {
        return Err(eyre::eyre!(
            "not a PPVM bytecode file (magic 0x{magic:08X})"
        ));
    }
    let version = read_u16(r)?;
    if version != PPVM_BYTECODE_VERSION {
        return Err(eyre::eyre!("unsupported bytecode version {version}"));
    }
    let header_size = read_u32(r)?;
    let n_qubits = read_u32(r)? as usize;
    let coefficient_threshold = read_f64(r)?;

    // Sections begin at `header_size`; skip any header bytes beyond v1's fixed
    // fields (forward compat / self-description).
    if header_size < HEADER_SIZE {
        return Err(eyre::eyre!(
            "header_size {header_size} smaller than minimum {HEADER_SIZE}"
        ));
    }
    skip_bytes(r, u64::from(header_size - HEADER_SIZE))?;

    // Don't pre-allocate from an untrusted count; grow as entries are read.
    let string_count = read_u32(r)?;
    let mut strings = Vec::new();
    for _ in 0..string_count {
        let len = read_u32(r)? as usize;
        let mut bytes = vec![0u8; len];
        r.read_exact(&mut bytes)?;
        strings.push(String::from_utf8(bytes)?);
    }

    let code_count = read_u32(r)?;
    let mut code = Vec::new();
    for _ in 0..code_count {
        code.push(PPVMInstruction::from_bytes(r)?);
    }

    Ok(PPVMModule {
        extra: PPVMDeviceInfo {
            magic,
            n_qubits,
            coefficient_threshold,
        },
        strings,
        code,
        ..Default::default()
    })
}

/// Serialize a module to an owned byte vector.
pub fn module_to_bytes(module: &PPVMModule) -> eyre::Result<Vec<u8>> {
    let mut buf = Vec::new();
    write_module(module, &mut buf)?;
    Ok(buf)
}

/// Reconstruct a module from a byte slice.
pub fn module_from_bytes(bytes: &[u8]) -> eyre::Result<PPVMModule> {
    read_module(&mut &bytes[..])
}

/// Cheap sniff: does this byte stream begin with the PPVM `.ssb` magic?
///
/// Reads the leading four bytes the same way [`read_module`] does — as a
/// little-endian `u32` — so a positive result here means [`read_module`] will
/// accept the magic. A stream shorter than the magic is not bytecode.
pub fn is_bytecode(bytes: &[u8]) -> bool {
    bytes.len() >= 4 && u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) == PPVM_MAGIC
}

/// "Dump": compile `.sst` source straight to the `.ssb` byte stream.
pub fn compile_to_bytes(source: &str) -> eyre::Result<Vec<u8>> {
    let module = crate::compile_program(source)?;
    module_to_bytes(&module)
}

fn read_u16<R: Read>(r: &mut R) -> eyre::Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}

fn read_u32<R: Read>(r: &mut R) -> eyre::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_f64<R: Read>(r: &mut R) -> eyre::Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

fn skip_bytes<R: Read>(r: &mut R, n: u64) -> eyre::Result<()> {
    let skipped = std::io::copy(&mut r.take(n), &mut std::io::sink())?;
    if skipped != n {
        return Err(eyre::eyre!("unexpected EOF skipping {n} header bytes"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_module() -> PPVMModule {
        PPVMModule::default()
    }

    #[test]
    fn round_trips_device_info() {
        let mut m = empty_module();
        m.extra.n_qubits = 7;
        m.extra.coefficient_threshold = 1e-9;

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn round_trips_code() {
        use crate::instruction::CircuitInstruction;
        use vihaco_cpu::Instruction as Cpu;

        let mut m = empty_module();
        m.extra.n_qubits = 2;
        m.code = vec![
            PPVMInstruction::Cpu(Cpu::Const(Value::U64(0))),
            PPVMInstruction::Circuit(CircuitInstruction::H),
            PPVMInstruction::Cpu(Cpu::Branch(1)),
            PPVMInstruction::Cpu(Cpu::ConditionalBranch(0, 1)),
            PPVMInstruction::Cpu(Cpu::Call(0, 1)),
            PPVMInstruction::Cpu(Cpu::Return(0)),
        ];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn read_honors_header_size_with_padding() {
        let mut m = empty_module();
        m.extra.n_qubits = 3;
        m.strings = vec!["hi".to_string()];
        m.code = vec![PPVMInstruction::Cpu(vihaco_cpu::Instruction::Return(0))];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();

        // Simulate a larger header: 4 padding bytes after the fixed fields,
        // with header_size bumped to match. The reader must skip to it.
        buf[6..10].copy_from_slice(&(HEADER_SIZE + 4).to_le_bytes());
        for i in 0..4 {
            buf.insert(HEADER_SIZE as usize + i, 0x00);
        }

        let back = read_module(&mut buf.as_slice()).unwrap();
        assert_eq!(back, m);
    }

    #[test]
    fn compile_to_bytes_round_trips_through_resolve() {
        let src = "device circuit.n_qubits 2;\n\
                   fn @main() {\n\
                       const.u64 0\n\
                       gate h\n\
                       const.u64 0\n\
                       const.u64 1\n\
                       gate cnot\n\
                       ret\n\
                   }\n";

        let bytes = compile_to_bytes(src).unwrap();
        let back = module_from_bytes(&bytes).unwrap();
        let expected = crate::compile_program(src).unwrap();

        assert_eq!(back, expected);
    }

    #[test]
    fn loaded_bytecode_executes_like_text() {
        let src = "device circuit.n_qubits 2;\n\
                   fn @main() {\n\
                       const.u64 0\n gate h\n\
                       const.u64 0\n const.u64 1\n gate cnot\n\
                       const.u64 0\n gate measure\n\
                       const.u64 1\n gate measure\n\
                       ret\n }\n";
        let bytes = compile_to_bytes(src).unwrap();

        let mut machine = crate::composite::PPVM::default();
        machine.load_bytecode(&bytes).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 2);
    }

    #[test]
    fn load_bytecode_file_reads_from_disk() {
        let src = "device circuit.n_qubits 1;\n\
                   fn @main() { const.u64 0\n gate measure\n ret }\n";
        let bytes = compile_to_bytes(src).unwrap();
        let path = std::env::temp_dir().join("ppvm_load_bytecode_file_test.ssb");
        std::fs::write(&path, &bytes).unwrap();

        let mut machine = crate::composite::PPVM::default();
        machine.load_bytecode_file(path.to_str().unwrap()).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_rejects_truncated_input() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        m.code = vec![PPVMInstruction::Cpu(vihaco_cpu::Instruction::Return(0))];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        buf.truncate(buf.len() - 3); // cut off mid-instruction

        assert!(read_module(&mut buf.as_slice()).is_err());
    }

    #[test]
    fn write_rejects_populated_functions_table() {
        use vihaco::module::{FunctionInfo, Signature};

        let mut m = empty_module();
        m.extra.n_qubits = 1;
        m.functions.push(FunctionInfo {
            name: 0,
            signature: Signature {
                params: vec![],
                ret: vec![],
            },
            local_count: 0,
            start_address: 0,
            end_address: 0,
            file: 0,
        });

        let mut buf = Vec::new();
        let err = write_module(&m, &mut buf).unwrap_err();
        assert!(err.to_string().contains("functions"), "err: {err}");
    }

    #[test]
    fn read_rejects_bad_magic() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        // Corrupt the magic (first 4 bytes).
        buf[0] ^= 0xFF;

        let err = read_module(&mut buf.as_slice()).unwrap_err();
        assert!(
            err.to_string().contains("not a PPVM bytecode file"),
            "err: {err}"
        );
    }

    #[test]
    fn round_trips_strings() {
        let mut m = empty_module();
        m.extra.n_qubits = 1;
        m.strings = vec![
            String::new(),
            "hello".to_string(),
            "tab\tnl\nquote\"".to_string(),
            "üñîçødé ⚛".to_string(),
        ];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn read_rejects_unsupported_version() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        // Bump the version (bytes 4..6) past the supported one.
        let bad = (PPVM_BYTECODE_VERSION + 1).to_le_bytes();
        buf[4] = bad[0];
        buf[5] = bad[1];

        let err = read_module(&mut buf.as_slice()).unwrap_err();
        assert!(
            err.to_string().contains("unsupported bytecode version"),
            "err: {err}"
        );
    }
}
