use clap::Parser;
use std::path::PathBuf;

use oboannotation::goannotation::process_file;

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
    let result = process_file(goa_path_str);
    match result {
        Ok(json_string) => println!("{}", json_string),
        Err(e) => eprint!("Could not extract GOA stats: {}", e.to_string())
    }

}