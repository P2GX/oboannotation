use clap::Parser;
use flate2::bufread::GzDecoder;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use oboannotation::{go::GoGafAnnotationLoader, io::AnnotationLoader};

#[derive(Parser, Debug)]
#[command(version = "0.1.6", about = "Oboannotation demo")]
struct Args {
    /// Path to the goa_human.gaf.gz file
    #[arg(long, value_name = "FILE")]
    goa: PathBuf,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let goa_path = args.goa;
    let goa_path_str = goa_path.to_str().expect("Path should include UTF-8 values");
    println!("processing {}", goa_path_str);

    let loader = GoGafAnnotationLoader;

    let read: Box<dyn Read> = if goa_path_str.ends_with(".gz") {
        // Decompress on the fly.
        Box::new(GzDecoder::new(BufReader::new(File::open(goa_path_str)?)))
    } else {
        Box::new(File::open(goa_path_str)?)
    };

    match loader.load_from_read(read) {
        Ok(annotations) => println!(
            "Loaded {:?} annotations. Skipped {:?} negated annotations",
            annotations.annotations.len(),
            annotations.negated_annotation_count
        ),
        Err(e) => eprintln!("Could not load GOA: {e}"),
    }

    Ok(())
}
