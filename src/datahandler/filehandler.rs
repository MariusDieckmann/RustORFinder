use std::{fs, path::Path};

pub fn read_fasta_normal(path: Box<Path>) -> String {
    fs::read_to_string(path).unwrap()
}
