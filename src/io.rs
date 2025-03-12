//! Load ontology annotations.
use std::fs::File;

use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Load annotations from path, [`Read`], or [`BufRead`].
pub trait AnnotationLoader<A> {
    /// Load annotation from a file path.
    fn load_from_path<P>(&self, path: P) -> anyhow::Result<A>
    where
        P: AsRef<Path>,
    {
        self.load_from_read(File::open(path)?)
    }

    /// Load annotation from a reader.
    fn load_from_read<R>(&self, read: R) -> anyhow::Result<A>
    where
        R: Read,
    {
        self.load_from_buf_read(BufReader::new(read))
    }

    /// Load annotation from a buffered reader.
    fn load_from_buf_read<R>(&self, read: R) -> anyhow::Result<A>
    where
        R: BufRead;
}
