use std::{io, sync::Mutex};

use super::outwriter::OutWriter;
use bio::io::gff;

pub struct GffWriter<T: io::Write> {
    writer: Mutex<gff::Writer<T>>,
}

impl<T: io::Write> GffWriter<T> {
    pub fn new(writer: T) -> Self {
        let gff_writer = gff::Writer::new(writer, gff::GffType::GFF3);
        return GffWriter {
            writer: Mutex::new(gff_writer),
        };
    }
}

impl<T: io::Write> OutWriter for GffWriter<T> {
    fn write(
        &self,
        orfs: crossbeam::channel::Receiver<crate::models::models::ORF>,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let mut gff_writer = self.writer.lock().unwrap();

        for orf in orfs {
            let orf_uuid = uuid::Uuid::new_v4();

            let mut gff_record = gff::Record::new();
            let start = gff_record.start_mut();
            *start = orf.start_position as u64;
            let stop = gff_record.end_mut();
            *stop = orf.stop_position as u64;
            let name = gff_record.seqname_mut();
            *name = orf_uuid.to_string();
            let strand = gff_record.strand_mut();
            *strand = match orf.direction {
                crate::models::models::Direction::FORWARD => "+".to_string(),
                crate::models::models::Direction::REVERSE => "-".to_string(),
            };
            let feature_type = gff_record.feature_type_mut();
            *feature_type = "ORF".to_string();

            gff_writer.write(&gff_record)?;
        }

        return Ok(());
    }
}
