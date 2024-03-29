use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    path::PathBuf,
};

use clap::Parser;
use datahandler::filehandler::SequenceFileType;
use rustyorffinder::{find_orfs, outwriter::outwriter::OutputType};

const MIN_NUM_THREADS: u8 = 4;

pub mod datahandler;
pub mod finder;
pub mod models;
pub mod outwriter;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// The path to the sequence file, can be any FD type, including e.g. /dev/stdin
    #[clap(short, long, value_name = "FILE")]
    sequence_path: Option<PathBuf>,

    /// Input kind
    #[clap(short = 't', long, arg_enum)]
    sequence_file_type: SequenceFileType,

    /// Number of threads
    #[clap(short, long, value_parser = validate_thread_number, default_value_t = 4)]
    num_threads: u8,

    /// Masked regions in GFF3 format
    #[clap(short, long, value_name = "FILE")]
    masked_gff3: Option<PathBuf>,

    /// Minimum size of an ORF
    #[clap(long, default_value_t = 30)]
    orf_min_length: usize,

    /// Output file format
    #[clap(short = 'h', long, arg_enum)]
    output_format: OutputType,

    /// Output file path
    #[clap(short = 'o', long)]
    output_file: Option<PathBuf>,
}

fn main() {
    env_logger::init();

    let cli = Cli::parse();
    let input_io: Box<dyn std::io::Read + Send + Sync> = match cli.sequence_path {
        Some(path) => {
            let file: Box<dyn std::io::Read + Send + Sync> =
                Box::new(OpenOptions::new().read(true).open(path).unwrap());
            file
        }
        None => Box::new(std::io::stdin()),
    };

    let sequences = datahandler::filehandler::read_input(input_io, cli.sequence_file_type).unwrap();

    let masked_areas = match cli.masked_gff3 {
        Some(value) => datahandler::mask_file::parse_mask_gff3_file(value).unwrap(),
        None => HashMap::new(),
    };

    let output_type = cli.output_format;

    let out_io: Box<dyn std::io::Write + Send + Sync> = match cli.output_file {
        Some(path) => Box::new(File::create(path).unwrap()),
        None => Box::new(std::io::stdout()),
    };

    let sequence = sequences.get(0).unwrap();
    find_orfs(
        sequence.to_string(),
        masked_areas,
        cli.num_threads,
        false,
        11,
        cli.orf_min_length,
        output_type,
        out_io,
    );
}

fn validate_thread_number(num_threads_string: &str) -> Result<u8, String> {
    let num_threads: u8 = num_threads_string
        .parse()
        .map_err(|_| format!("`{}` isn't a valid number of threads", num_threads_string))?;

    if num_threads < MIN_NUM_THREADS {
        return Err(format!(
            "Insufficient number of threads provided: At least {} threads required to run.",
            MIN_NUM_THREADS
        ));
    }

    return Ok(num_threads);
}
