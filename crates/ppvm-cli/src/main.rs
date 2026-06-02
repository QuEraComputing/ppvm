use clap::{Parser, Subcommand};
use eyre::Result;

mod commands;

#[derive(Parser)]
#[command(name = "ppvm")]
#[command(about = "Pauli propagation virtual machine", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Parse a .sst file and output the AST
    Parse {
        /// Input .sst file
        #[arg(value_name = "FILE")]
        file: String,

        /// Output format
        #[arg(short, long, value_enum, default_value = "pretty")]
        format: commands::Format,
    },
    /// Compile a .sst file to bytecode
    Dump {
        /// Input .sst file
        #[arg(value_name = "FILE")]
        file: String,

        /// Output file (optional, defaults to <file_name>.ssb)
        #[arg(short, long)]
        output: Option<String>,

        /// Overwrite the output file if it already exists
        #[arg(short, long)]
        force: bool,
    },

    /// Run a .sst or .ssb program
    Run {
        /// Input file (.sst source or .ssb bytecode)
        #[arg(value_name = "FILE")]
        file: String,

        /// Suppress the measurement record
        #[arg(short, long)]
        quiet: bool,

        /// Measurement output format
        #[arg(short, long, value_enum, default_value = "bits")]
        format: commands::MeasurementFormat,
    },

    /// Step through a program interactively, pausing at `breakpoint` instructions
    Debug {
        /// Input file (.sst source or .ssb bytecode)
        #[arg(value_name = "FILE")]
        file: String,

        /// Pause before the first instruction so any program can be stepped
        #[arg(short, long)]
        break_at_start: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Parse { file, format } => {
            commands::parse(&file, format)?;
        }
        Commands::Dump {
            file,
            output,
            force,
        } => {
            commands::dump(&file, output.as_deref(), force)?;
        }
        Commands::Run {
            file,
            quiet,
            format,
        } => {
            commands::run(&file, quiet, format)?;
        }
        Commands::Debug {
            file,
            break_at_start,
        } => {
            commands::debug(&file, break_at_start)?;
        }
    }

    Ok(())
}
