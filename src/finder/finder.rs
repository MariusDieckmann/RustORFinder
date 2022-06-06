use std::result::Result;

use crossbeam::channel::Receiver;

use super::threaded_finder::ORF;

trait Finder {}

trait ORFWriter {
    fn write_orfs(orfs_chan_recv: Receiver<ORF>) -> Result<()>;
}
