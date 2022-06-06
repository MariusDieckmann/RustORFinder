use std::io::Write;

use crossbeam::{
    channel::{bounded, Receiver, Sender},
    thread,
};
use serde::{Deserialize, Serialize};

pub struct ThreadedFinder {
    pub fw_sequence: String,
    pub rev_sequence: String,
    pub start_codons: Vec<String>,
    pub stop_codons: Vec<String>,
    pub min_len: usize,
}

pub struct ORFPositions {
    pub start_positions: Vec<usize>,
    pub stop_position: usize,
    pub strand: Direction,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct ORF {
    start_position: usize,
    stop_position: usize,
    sequence: String,
    direction: Direction,
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub enum Direction {
    FORWARD,
    REVERSE,
}

impl ThreadedFinder {
    pub fn new(sequence: String) -> Self {
        let rev_seq: String = sequence
            .chars()
            .rev()
            .map(|x| match x {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => panic!("unexpected value found: {}", x),
            })
            .collect();
        let finder = ThreadedFinder {
            fw_sequence: sequence,
            rev_sequence: rev_seq,
            start_codons: vec!["ATG".to_string()],
            stop_codons: vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()],
            min_len: 30,
        };

        return finder;
    }

    pub fn find_orfs(&self) {
        thread::scope(|s| {
            let (sender, recv) = bounded(1000);
            let (transcribe_sender, transcribe_recv) = bounded(1000);

            let mut task_handles = Vec::new();
            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(0, Direction::FORWARD, cloned_sender));
            task_handles.push(handle);

            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(1, Direction::FORWARD, cloned_sender));
            task_handles.push(handle);

            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(2, Direction::FORWARD, cloned_sender));
            task_handles.push(handle);

            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(0, Direction::REVERSE, cloned_sender));
            task_handles.push(handle);

            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(1, Direction::REVERSE, cloned_sender));
            task_handles.push(handle);

            let cloned_sender = sender.clone();
            let handle = s.spawn(|_| self.orf_reader(2, Direction::REVERSE, cloned_sender));
            task_handles.push(handle);

            for _ in 0..10 {
                let cloned_recv_0 = recv.clone();
                let cloned_transcribe_sender = transcribe_sender.clone();
                let transcriber =
                    s.spawn(|_| self.transcribe_orfs(cloned_recv_0, cloned_transcribe_sender));
                task_handles.push(transcriber);
            }

            let cloned_transcribe_recv = transcribe_recv.clone();
            let writer = s.spawn(|_| self.write_out(cloned_transcribe_recv));
            task_handles.push(writer);

            drop(sender);
            drop(recv);
            drop(transcribe_sender);
            drop(transcribe_recv);

            for handle in task_handles {
                handle.join().unwrap();
            }
        })
        .unwrap();
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
            if self.start_codons.contains(&codon) {
                starts.push(n as usize);
            }

            if self.stop_codons.contains(&codon) {
                let orf = ORFPositions {
                    start_positions: starts,
                    stop_position: n as usize,
                    strand: direction,
                };

                starts = Vec::new();

                sender.send(orf).unwrap();
            }
        }
    }

    fn transcribe_orfs(&self, recv: Receiver<ORFPositions>, sender: Sender<Vec<u8>>) {
        for orf_positions in recv.iter() {
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

                let mut orf = serde_json::to_string(&orf).unwrap().as_bytes().to_vec();
                orf.push(10);
                sender.send(orf).unwrap();
            }
        }
    }

    fn write_out(&self, recv: Receiver<Vec<u8>>) {
        let mut i = 0;
        for orf in recv.iter() {
            i = i + 1;
            //std::io::stdout().write_all(orf.as_slice()).unwrap();
        }
        println!("{}", i);
    }
}

#[cfg(test)]
mod tests {
    use super::ThreadedFinder;

    #[test]
    fn test_orf_reader() {
        let sequence = "ATG";
        let finder = ThreadedFinder::new(sequence.to_string());
        finder.find_orfs();
    }
}
