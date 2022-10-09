use bio::io::fasta;
use clap::ValueEnum;
use std::{fs, path::PathBuf};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum SequenceFileType {
    RawSequence,
    Fasta,
}

pub fn read_input(
    path: Box<PathBuf>,
    file_type: SequenceFileType,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let sequences = match file_type {
        SequenceFileType::RawSequence => read_fasta_normal(path)?,
        SequenceFileType::Fasta => read_fasta(path)?,
    };
    return Ok(sequences);
}

fn read_fasta_normal(path: Box<PathBuf>) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let sequence = fs::read_to_string(path.as_path())?;

    return Ok(vec![sequence]);
}

fn read_fasta(path: Box<PathBuf>) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let fasta_file = std::fs::File::open(path.as_path()).unwrap();
    let fasta_reader = fasta::Reader::new(fasta_file);

    let mut sequences = Vec::new();
    for record in fasta_reader.records() {
        let fasta_entry = record?.seq().to_vec();
        let sequence = String::from_utf8(fasta_entry)?;
        sequences.push(sequence);
    }

    return Ok(sequences);
}
