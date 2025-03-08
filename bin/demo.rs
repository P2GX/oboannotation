use clap::Parser;
use std::path::PathBuf;

use oboannotation::goannotation::GoAnnotations;

#[derive(Parser, Debug)]
#[command(version = "0.1.6", about = "Fenominal implementation in Rust")]
struct Args {
   /// Path to the goa_human.gaf.gz file
   #[arg(long, value_name = "FILE")]
   goa: PathBuf,

}



fn main() {
    let args = Args::parse();
    let goa_path = args.goa;
    let goa_path_str: &str = goa_path.to_str().expect("Invalid UTF-8 in path");
    println!("processing {}", goa_path_str);
    let goannots = GoAnnotations::new(goa_path_str).expect("Could not create GO Annotations");
    let result = goannots.get_annotation_statistics_json();
    match result {
        Ok(json_string) => println!("{}", json_string),
        Err(e) => eprint!("Could not extract GOA stats: {}", e.to_string())
    }
    // now count the annotations for the first ten genes
    println!("##############\nShowing the first ten GO annotation profiles");
    let annot_map = goannots.get_annotation_map();
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        let count = term_ids.len();
        println!("{} has {} unique GO annotations", symbol, count);
        if c > 10 {
            break;
        } else {
            c += 1;
        }
    }

}