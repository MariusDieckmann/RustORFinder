use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct ORF {
    pub start_position: usize,
    pub stop_position: usize,
    pub sequence: String,
    pub direction: Direction,
}

pub struct ORFPositions {
    pub start_positions: Vec<usize>,
    pub stop_position: usize,
    pub strand: Direction,
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub enum Direction {
    FORWARD,
    REVERSE,
}
