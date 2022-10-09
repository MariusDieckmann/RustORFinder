use std::collections::HashMap;

use crate::{
    datahandler::{self, trans_table::TranslationalTable},
    models::models::{Direction, ORFPositions, ORF},
    outwriter::outwriter::OutWriter,
};
use crossbeam::channel::{bounded, Receiver, Sender};
use crossbeam::thread;

enum NextThreadType {
    Finder,
    Transcriber,
}

pub struct ThreadedFinder {
    pub fw_sequence: String,
    pub rev_sequence: String,
    pub masked_areas: HashMap<u64, u64>,
    pub translational_table: TranslationalTable,
    pub outwriter: Box<dyn OutWriter + Send + Sync + 'static>,
    pub min_len: usize,
    pub circular: bool,
}

impl ThreadedFinder {
    pub fn new(
        sequence: String,
        translational_table_id: u8,
        masked_areas: HashMap<u64, u64>,
        circular: bool,
        outwriter: Box<dyn OutWriter + Send + Sync + 'static>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let table =
            datahandler::trans_table::parse_translational_table(translational_table_id)?.unwrap();

        let rev_seq: String = sequence.chars().rev().collect();
        let finder = ThreadedFinder {
            fw_sequence: sequence,
            rev_sequence: rev_seq,
            translational_table: table,
            masked_areas: masked_areas,
            min_len: 30,
            circular: circular,
            outwriter: outwriter,
        };

        return Ok(finder);
    }

    pub fn run(&self, num_threads: u8) {
        let (orf_positions_sender, orf_positions_recv) = bounded(100);
        let (orf_sender, orf_receiver) = bounded(100);
        let (orf_reader_start_pos_sender, orf_reader_start_pos_recv) = bounded(6);

        for i in 0..3 {
            orf_reader_start_pos_sender
                .send((i, Direction::FORWARD))
                .unwrap();
            orf_reader_start_pos_sender
                .send((i, Direction::REVERSE))
                .unwrap();
        }

        let mut next_thread_type = NextThreadType::Finder;

        thread::scope(|s| {
            let cloned_orf_receiver = orf_receiver.clone();
            s.spawn(|_| self.outwriter.write(cloned_orf_receiver));
            for _ in 0..num_threads - 1 {
                match next_thread_type {
                    NextThreadType::Finder => {
                        let cloned_orf_reader_start_pos_recv = orf_reader_start_pos_recv.clone();
                        let cloned_orf_positions_sender = orf_positions_sender.clone();

                        s.spawn(|_| {
                            self.find_orfs(
                                cloned_orf_reader_start_pos_recv,
                                cloned_orf_positions_sender,
                            )
                        });
                        next_thread_type = NextThreadType::Transcriber;
                    }
                    NextThreadType::Transcriber => {
                        let cloned_orf_positions_recv = orf_positions_recv.clone();
                        let cloned_orf_sender = orf_sender.clone();

                        s.spawn(|_| {
                            self.transcribe_orfs(cloned_orf_positions_recv, cloned_orf_sender)
                        });

                        next_thread_type = NextThreadType::Finder;
                    }
                }
            }

            drop(orf_reader_start_pos_sender);
            drop(orf_reader_start_pos_recv);
            drop(orf_positions_sender);
            drop(orf_positions_recv);
            drop(orf_sender);
            drop(orf_receiver);
        })
        .unwrap();
    }

    fn find_orfs(
        &self,
        orf_reader_start_pos_recv: Receiver<(u64, Direction)>,
        orf_positions_sender: Sender<ORFPositions>,
    ) {
        for start_pos in orf_reader_start_pos_recv {
            self.orf_reader(start_pos.0, start_pos.1, orf_positions_sender.clone());
        }
    }

    fn orf_reader(&self, offset: u64, direction: Direction, sender: Sender<ORFPositions>) {
        let sequence: &String;
        match direction {
            Direction::FORWARD => sequence = &self.fw_sequence,
            Direction::REVERSE => sequence = &self.rev_sequence,
        }

        let mut starts = Vec::new();

        for n in (offset..sequence.len() as u64).step_by(3) {
            if n + 2 >= sequence.len() as u64 {
                break;
            }

            let i = n as usize;
            let codon = sequence[i..i + 3].to_string();
            if self.translational_table.start_codons.contains(&codon) {
                starts.push(n as usize);
            }

            if self.translational_table.stop_codons.contains(&codon) {
                let orf = ORFPositions {
                    start_positions: starts,
                    stop_position: n as usize,
                    strand: direction,
                };

                starts = Vec::new();

                sender.send(orf).unwrap();
            }
        }

        if self.circular {
            for n in (offset..sequence.len() as u64).step_by(3) {
                if n + 2 >= sequence.len() as u64 {
                    break;
                }

                let i = n as usize;
                let codon = sequence[i..i + 3].to_string();
                if self.translational_table.stop_codons.contains(&codon) {
                    let orf = ORFPositions {
                        start_positions: starts,
                        stop_position: n as usize,
                        strand: direction,
                    };

                    sender.send(orf).unwrap();

                    break;
                }
            }
        }
    }

    fn transcribe_orfs(&self, orf_positions_recv: Receiver<ORFPositions>, orf_sender: Sender<ORF>) {
        for orf_positions in orf_positions_recv.iter() {
            let sequence: String;
            match orf_positions.strand {
                Direction::FORWARD => sequence = self.fw_sequence.clone(),
                Direction::REVERSE => sequence = self.rev_sequence.clone(),
            }

            for start_pos in orf_positions.start_positions {
                let diff = orf_positions.stop_position - start_pos;
                if diff <= self.min_len {
                    continue;
                }

                let subsequence = &sequence.as_bytes()[start_pos..orf_positions.stop_position];
                let subsequence_string = std::str::from_utf8(subsequence).unwrap().to_string();

                let orf = ORF {
                    start_position: start_pos as usize,
                    stop_position: orf_positions.stop_position as usize,
                    sequence: subsequence_string,
                    direction: orf_positions.strand,
                };

                orf_sender.send(orf).unwrap();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs};

    use crate::outwriter::count_writer::CountWriter;

    use super::ThreadedFinder;

    #[test]
    fn threaded_finder_full() {
        let count_writer = Box::new(CountWriter {});
        let masked_areas = HashMap::new();
        let sequence = fs::read_to_string("resources/sequences/NC_011604.1_normal.fasta").unwrap();
        let threaded_finder =
            ThreadedFinder::new(sequence, 11, masked_areas, false, count_writer).unwrap();
        threaded_finder.run(4);
    }
}
