use serde::{Deserialize, Serialize};

const TRANS_TABLE_11: &'static [u8] =
    include_bytes!("../../resources/codon_tables/trans_table_11.json");

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TranslationalTable {
    pub start_codons: Vec<String>,
    pub stop_codons: Vec<String>,
}

pub fn parse_translational_table(
    trans_table_number: u8,
) -> Result<Option<TranslationalTable>, Box<dyn std::error::Error>> {
    let codon_translation_table: TranslationalTable = match trans_table_number {
        11 => serde_json::from_slice(TRANS_TABLE_11)?,
        _ => return Ok(None),
    };

    return Ok(Some(codon_translation_table));
}

#[cfg(test)]
mod tests {
    use super::parse_translational_table;

    #[test]
    fn test_parse_translational_table_11() {
        parse_translational_table(11).unwrap().unwrap();
    }
}
