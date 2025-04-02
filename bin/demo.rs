use anyhow::bail;
use clap::Parser;
use flate2::read::GzDecoder;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use oboannotation::{
    go::{GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};

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

    println!("processing {:?}", goa_path);

    let reader: Box<dyn BufRead> = open_for_reading(goa_path)?;
    
    let loader = GoGafAnnotationLoader;
    let annotations = match loader.load_from_buf_read(reader) {
        Ok(annotations) => {
            println!(
                "Loaded {:?} annotations. Skipped {:?} negated annotations",
                annotations.annotations.len(),
                annotations.negated_annotation_count
            );
            annotations
        }
        Err(e) => {
            bail!("Could not load GOA: {e}")
        }
    };
    // now count the annotations for the first ten genes
    println!("##############\nShowing the first ten GO annotation profiles");
    let annot_map = get_annotation_map(&annotations);
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        println!("{} has {} unique GO annotations", symbol, term_ids.len());
        if c > 10 {
            break;
        } else {
            c += 1;
        }
    }

    Ok(())
}

fn open_for_reading<P: AsRef<Path>>(goa_path: P) -> anyhow::Result<Box<dyn BufRead>> {
    Ok(if let Some(extension) = goa_path.as_ref().extension() {
        if extension == "gz" {
            Box::new(BufReader::new(GzDecoder::new(File::open(goa_path)?)))
        } else {
            Box::new(BufReader::new(File::open(goa_path)?))
        }
    } else {
        Box::new(BufReader::new(File::open(goa_path)?))
    })
}
