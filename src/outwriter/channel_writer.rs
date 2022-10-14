use std::{
    cell::RefCell,
    sync::{Mutex, RwLock},
    thread::{self, Thread},
};

use crossbeam::channel::{Receiver, Sender};

use crate::models::models::ORF;

use super::outwriter::OutWriter;

pub struct ChannelWriter {
    send: Sender<ORF>,
    recv: Receiver<ORF>,
}

impl ChannelWriter {
    pub fn new(send: Sender<ORF>, recv: Receiver<ORF>) -> Self {
        return ChannelWriter {
            send: send,
            recv: recv,
        };
    }
}

impl OutWriter for ChannelWriter {
    fn write(&self, orfs: Receiver<ORF>) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        for orf in orfs {
            self.send.send(orf).unwrap();
        }

        Ok(())
    }
}
