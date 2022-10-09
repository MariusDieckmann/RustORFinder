use std::io::Write;

use crossbeam::channel::Receiver;

use crate::models::models::ORF;
use clap::ValueEnum;

use super::{count_writer::CountWriter, gff_writer::GffWriter};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum OutputType {
    Counter,
    GFF3,
}

pub trait OutWriter {
    fn write(&self, orfs: Receiver<ORF>) -> Result<(), Box<dyn std::error::Error + Send + Sync>>;
}

pub fn get_writer(
    output_format: OutputType,
    out: Box<dyn Write + Send + Sync>,
) -> Box<dyn OutWriter + Send + Sync + 'static> {
    let writer: Box<dyn OutWriter + Send + Sync> = match output_format {
        OutputType::Counter => Box::new(CountWriter {}),
        OutputType::GFF3 => Box::new(GffWriter::new(out)),
    };

    return writer;
}
