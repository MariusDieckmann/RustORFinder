use std::path::Path;

use datahandler::filehandler;

pub mod datahandler;
pub mod finder;

fn main() {
    let data = filehandler::read_fasta_normal(
        Box::new(
            Path::new("/Users/mariusdieckmann/Data/orffinder-ng/NC_002695.2_normal.fna").to_owned(),
        )
        .into_boxed_path(),
    );

    let finder = finder::threaded_finder::ThreadedFinder::new(data);
    //println!("{:?}", finder.rev_sequence)
    finder.find_orfs();
}
