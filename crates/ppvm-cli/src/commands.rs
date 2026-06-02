use eyre::Result;
use ppvm_vihaco::run_file;
use std::path::Path;

pub fn run(file: &str, show_measurements: bool) -> Result<()> {
    let ppvm = run_file(file)?;
    if show_measurements {
        let outcomes = ppvm.measurement_record();
        println!("Measurement record:\n{:?}", outcomes);
    } else {
        println!("Successfully ran file {}", file);
    }
    Ok(())
}

pub fn parse(file: &str, format: &str) -> Result<()> {
    let source = std::fs::read_to_string(file)?;
    let parsed = ppvm_vihaco::parse_program(&source)?;

    match format {
        "json" => {
            eprintln!("Warning: JSON format not yet supported for AST, using debug format");
            println!("{:#?}", parsed);
        }
        "debug" => {
            println!("{:#?}", parsed);
        }
        _ => {
            // Pretty summary.
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

pub fn dump(file: &str, output: Option<&str>) -> Result<()> {
    let output_file = match output {
        Some(output_file_name) => output_file_name.to_string(),
        None => Path::new(file)
            .with_extension("ssb")
            .to_string_lossy()
            .into_owned(),
    };

    ppvm_vihaco::dump_file(file, &output_file)?;
    eprintln!("Bytecode written to {output_file}");
    Ok(())
}
