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
