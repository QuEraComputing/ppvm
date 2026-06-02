use eyre::{Result, WrapErr};
use ppvm_vihaco::measurements::{MeasurementOutcome, MeasurementResult};
use ppvm_vihaco::run_file;
use std::path::Path;

/// Output format for `parse`.
#[derive(Clone, Debug, clap::ValueEnum)]
pub enum Format {
    Pretty,
    Debug,
    Json,
}

/// Output format for the measurement record from `run`.
#[derive(Clone, Debug, clap::ValueEnum)]
pub enum MeasurementFormat {
    /// One bit string per measurement event, space-separated; lost qubit = `L`.
    Bits,
    /// Raw debug representation of the record.
    Debug,
}

pub fn run(file: &str, quiet: bool, format: MeasurementFormat) -> Result<()> {
    let ppvm = run_file(file).wrap_err_with(|| format!("failed to run {file}"))?;
    if quiet {
        return Ok(());
    }
    let record = ppvm.measurement_record();
    match format {
        MeasurementFormat::Bits => println!("Measurements: {}", format_bits(&record)),
        MeasurementFormat::Debug => println!("Measurement record:\n{:?}", record),
    }
    Ok(())
}

/// Render each measurement event as a bit string (lost qubit = `L`), events
/// space-separated. Empty record renders as `(none)`.
fn format_bits(record: &[MeasurementResult]) -> String {
    if record.is_empty() {
        return "(none)".to_string();
    }
    record
        .iter()
        .map(|event| {
            event
                .iter()
                .map(|outcome| match outcome {
                    MeasurementOutcome::Zero => '0',
                    MeasurementOutcome::One => '1',
                    MeasurementOutcome::Lost => 'L',
                })
                .collect::<String>()
        })
        .collect::<Vec<_>>()
        .join(" ")
}

pub fn parse(file: &str, format: Format) -> Result<()> {
    let source =
        std::fs::read_to_string(file).wrap_err_with(|| format!("failed to read {file}"))?;
    let parsed = ppvm_vihaco::parse_program(&source)?;

    match format {
        Format::Json => {
            eprintln!("Warning: JSON format not yet supported for AST, using debug format");
            println!("{:#?}", parsed);
        }
        Format::Debug => {
            println!("{:#?}", parsed);
        }
        Format::Pretty => {
            println!("Module:");
            println!("  Headers: {}", parsed.headers.len());
            for (i, header) in parsed.headers.iter().enumerate() {
                println!("    [{}] {:?}", i, header);
            }
            println!("  Functions: {}", parsed.functions.len());
            for (i, func) in parsed.functions.iter().enumerate() {
                println!(
                    "    [{}] {}({} params, {} body items)",
                    i,
                    func.name,
                    func.params.len(),
                    func.body.len()
                );
            }
        }
    }

    Ok(())
}

pub fn dump(file: &str, output: Option<&str>, force: bool) -> Result<()> {
    let output_file = match output {
        Some(output_file_name) => output_file_name.to_string(),
        None => Path::new(file)
            .with_extension("ssb")
            .to_string_lossy()
            .into_owned(),
    };

    // Don't clobber an existing file unless asked to.
    if !force && Path::new(&output_file).exists() {
        return Err(eyre::eyre!(
            "{output_file} already exists; pass --force to overwrite"
        ));
    }

    ppvm_vihaco::dump_file(file, &output_file)
        .wrap_err_with(|| format!("failed to dump {file}"))?;
    eprintln!("Bytecode written to {output_file}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Minimal program that compiles and measures q0 in |0> (deterministic).
    const PROGRAM: &str =
        "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n gate measure\n ret }\n";

    fn row(outcomes: &[MeasurementOutcome]) -> MeasurementResult {
        outcomes.iter().copied().collect()
    }

    /// Write `contents` to a uniquely-named temp file and return its path.
    fn temp_file(name: &str, contents: &str) -> String {
        let path = std::env::temp_dir().join(name);
        fs::write(&path, contents).unwrap();
        path.to_string_lossy().into_owned()
    }

    // ─── format_bits ───────────────────────────────────────────────────

    #[test]
    fn format_bits_empty_record_is_none() {
        assert_eq!(format_bits(&[]), "(none)");
    }

    #[test]
    fn format_bits_concatenates_qubits_within_an_event() {
        let record = vec![row(&[
            MeasurementOutcome::One,
            MeasurementOutcome::Zero,
            MeasurementOutcome::One,
        ])];
        assert_eq!(format_bits(&record), "101");
    }

    #[test]
    fn format_bits_separates_events_with_spaces() {
        let record = vec![
            row(&[MeasurementOutcome::One]),
            row(&[MeasurementOutcome::Zero]),
        ];
        assert_eq!(format_bits(&record), "1 0");
    }

    #[test]
    fn format_bits_renders_lost_qubit_as_l() {
        let record = vec![
            row(&[MeasurementOutcome::One, MeasurementOutcome::Lost]),
            row(&[MeasurementOutcome::Zero]),
        ];
        assert_eq!(format_bits(&record), "1L 0");
    }

    // ─── run ───────────────────────────────────────────────────────────

    #[test]
    fn run_succeeds_on_valid_file() {
        let src = temp_file("ppvm_cli_run_ok.sst", PROGRAM);
        let res = run(&src, true, MeasurementFormat::Bits);
        let _ = fs::remove_file(&src);
        assert!(res.is_ok(), "got: {res:?}");
    }

    #[test]
    fn run_errors_with_context_on_missing_file() {
        let err = run("/no/such/file.sst", false, MeasurementFormat::Bits).unwrap_err();
        assert!(err.to_string().contains("failed to run"), "got: {err}");
    }

    // ─── parse ─────────────────────────────────────────────────────────

    #[test]
    fn parse_succeeds_on_valid_file() {
        let src = temp_file("ppvm_cli_parse_ok.sst", PROGRAM);
        let res = parse(&src, Format::Debug);
        let _ = fs::remove_file(&src);
        assert!(res.is_ok(), "got: {res:?}");
    }

    #[test]
    fn parse_errors_with_context_on_missing_file() {
        let err = parse("/no/such/file.sst", Format::Pretty).unwrap_err();
        assert!(err.to_string().contains("failed to read"), "got: {err}");
    }

    // ─── dump ──────────────────────────────────────────────────────────

    #[test]
    fn dump_writes_default_ssb_path_when_output_omitted() {
        let src = temp_file("ppvm_cli_dump_default.sst", PROGRAM);
        let expected = Path::new(&src).with_extension("ssb");
        let _ = fs::remove_file(&expected); // clear any leftover from a prior run

        dump(&src, None, false).unwrap();
        assert!(expected.exists(), "default .ssb should have been written");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&expected);
    }

    #[test]
    fn dump_writes_explicit_output_path() {
        let src = temp_file("ppvm_cli_dump_explicit.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_explicit_out.ssb");
        let _ = fs::remove_file(&out);

        dump(&src, Some(out.to_str().unwrap()), false).unwrap();
        assert!(out.exists());

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    #[test]
    fn dump_refuses_to_clobber_without_force() {
        let src = temp_file("ppvm_cli_dump_clobber.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_clobber_out.ssb");
        fs::write(&out, b"existing").unwrap();

        let err = dump(&src, Some(out.to_str().unwrap()), false).unwrap_err();
        assert!(err.to_string().contains("already exists"), "got: {err}");
        // The existing file must be left untouched.
        assert_eq!(fs::read(&out).unwrap(), b"existing");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    #[test]
    fn dump_overwrites_existing_with_force() {
        let src = temp_file("ppvm_cli_dump_force.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_force_out.ssb");
        fs::write(&out, b"existing").unwrap();

        dump(&src, Some(out.to_str().unwrap()), true).unwrap();
        // Replaced with real bytecode, not the placeholder.
        assert_ne!(fs::read(&out).unwrap(), b"existing");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }
}
