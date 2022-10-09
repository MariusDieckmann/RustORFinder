use super::outwriter::OutWriter;

pub struct CountWriter {}

impl OutWriter for CountWriter {
    fn write(
        &self,
        orfs: crossbeam::channel::Receiver<crate::models::models::ORF>,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let mut i: u64 = 0;
        for _ in orfs {
            i = i + 1;
        }

        println!("Count: {}", i);
        return Ok(());
    }
}
