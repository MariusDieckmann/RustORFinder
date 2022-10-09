use std::{collections::HashMap, path::PathBuf};

use bio::io::gff;

pub fn parse_mask_gff3_file(
    path: PathBuf,
) -> Result<HashMap<u64, u64>, Box<dyn std::error::Error>> {
    let gff3_file = std::fs::File::open(path.as_path())?;
    let mut gff_reader = gff::Reader::new(gff3_file, gff::GffType::GFF3);
    let mut mask_starts = HashMap::new();

    for record_result in gff_reader.records() {
        let record = record_result?;
        mask_starts.insert(record.start().to_owned(), record.end().to_owned());
    }

    return Ok(mask_starts);
}
