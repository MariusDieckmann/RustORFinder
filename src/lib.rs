use std::collections::HashMap;

use outwriter::outwriter::{get_writer, OutputType};

pub mod datahandler;
pub mod finder;
mod models;
pub mod outwriter;

/// Finds all open reading frame (ORFs) using a multithreaded approach.
/// At least 4 threads are required.
/// * `sequence` - The sequence on which to search ORFs on
/// * `masked_areas` - A map of masked areas in the sequence that is not relevant
/// * `threads` - Number of threads to use, at least 4 are required
/// * `circular` - Indicates if the sequence is circular
/// * `trans_table` - Number of the translational table to use - this will change in a future update
/// * `out_format` - The output format in which the results should be written in
/// * `out_target`- The output target which to write the results to
pub fn find_orfs(
    sequence: String,
    masked_areas: HashMap<u64, u64>,
    threads: u8,
    circular: bool,
    trans_table: u8,
    min_len: usize,
    out_format: OutputType,
    out_target: Box<dyn std::io::Write + Send + Sync>,
) {
    let writer = get_writer(out_format, out_target);

    let finder = finder::threaded_finder::ThreadedFinder::new(
        sequence.to_string(),
        trans_table,
        masked_areas.clone(),
        circular,
        min_len,
        writer,
    )
    .unwrap();
    finder.run(threads);
}
