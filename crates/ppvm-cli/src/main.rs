// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use clap::{Parser, Subcommand};
use eyre::Result;

mod commands;

#[derive(Parser)]
#[command(name = "ppvm")]
#[command(about = "Pauli propagation virtual machine", long_about = None)]
pub struct Cli {
    /// Number of threads for all parallel work (1 = fully serial & deterministic)
    #[arg(short, long, default_value_t = 1, value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..))]
    threads: usize,

    /// Subcommand to run.
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

        /// Number of shots to run
        #[arg(short, long, default_value = "1")]
        shots: usize,

        /// Seed the RNG for reproducible results
        #[arg(long)]
        seed: Option<u64>,

        /// Write results to a file instead of stdout
        #[arg(short, long)]
        output: Option<String>,

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

    // Size the process-wide thread pool once; governs all parallelism (across
    // shots and within a single machine). `--threads 1` is fully serial.
    ppvm_vihaco::shots::set_global_threads(cli.threads)?;

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
            shots,
            seed,
            output,
            quiet,
            format,
        } => {
            commands::run(&file, shots, seed, output.as_deref(), quiet, format)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_zero_threads() {
        // `--threads 0` must be rejected (0 would break pool sizing). The
        // subcommand is valid, so a zero thread count is the only parse error.
        assert!(Cli::try_parse_from(["ppvm", "--threads", "0", "parse", "f.sst"]).is_err());
    }

    #[test]
    fn accepts_positive_and_default_threads() {
        assert!(Cli::try_parse_from(["ppvm", "--threads", "1", "parse", "f.sst"]).is_ok());
        assert!(Cli::try_parse_from(["ppvm", "parse", "f.sst"]).is_ok());
    }
}
