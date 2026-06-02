use clap::{Parser, Subcommand};
use eyre::Result;
use ppvm_cli::commands;

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

        /// Output format (json, debug, pretty)
        #[arg(short, long, default_value = "pretty")]
        format: String,
    },
    /// Compile a .sst file to bytecode
    Dump {
        /// Input .sst file
        #[arg(value_name = "FILE")]
        file: String,

        /// Output file (optional, defaults to <file_name>.ssb)
        #[arg(short, long)]
        output: Option<String>,
    },

    /// Run a .sst or .ssb program
    Run {
        /// Input file (.sst source or .ssb bytecode)
        #[arg(value_name = "FILE")]
        file: String,

        /// Show measurement output
        #[arg(short, long, default_value = "true")]
        show_measurements: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Parse { file, format } => {
            commands::parse(&file, &format)?;
        }
        Commands::Dump { file, output } => {
            commands::dump(&file, output.as_deref())?;
        }
        Commands::Run {
            file,
            show_measurements,
        } => {
            commands::run(&file, show_measurements)?;
        }
    }

    Ok(())
}
