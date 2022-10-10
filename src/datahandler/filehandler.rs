use bio::io::fasta;
use clap::ValueEnum;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum SequenceFileType {
    RawSequence,
    Fasta,
}

pub fn read_input(
    input: Box<dyn std::io::Read + Send + Sync>,
    file_type: SequenceFileType,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let sequences = match file_type {
        SequenceFileType::RawSequence => read_fasta_normal(input)?,
        SequenceFileType::Fasta => read_fasta(input)?,
    };
    return Ok(sequences);
}

fn read_fasta_normal(
    mut input: Box<dyn std::io::Read + Send + Sync>,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut string = String::new();
    input.read_to_string(&mut string)?;

    return Ok(vec![string]);
}

fn read_fasta(
    input: Box<dyn std::io::Read + Send + Sync>,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let fasta_reader = fasta::Reader::new(input);

    let mut sequences = Vec::new();
    for record in fasta_reader.records() {
        let fasta_entry = record?.seq().to_vec();
        let sequence = String::from_utf8(fasta_entry)?;
        sequences.push(sequence);
    }

    return Ok(sequences);
}
