//! Load ontology annotations.
use std::fs::File;

use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationIssue {
    /// Line number.
    pub n: usize,
    /// The issue description
    pub reason: String,
}

impl ValidationIssue {
    /// Create a new validation issue for the `n`th line.
    pub fn new<T>(n: usize, reason: &T) -> Self
    where
        T: ToString,
    {
        Self {
            n,
            reason: reason.to_string(),
        }
    }
}

/// Error with issues encountered in annotation parsing.
#[derive(Error, Debug)]
pub enum AnnotationLoadError {
    /// Error related to I/O from the underlying data source.
    #[error("I/O error")]
    IoError(#[from] io::Error),
    /// Error reporting failure of data validation (e.g. a misformatted CURIE string).
    #[error("Validation issues encountered")]
    ValidationError(Vec<ValidationIssue>),
    #[error("Error: {0}")]
    Error(String),
}

/// Load annotations from path, [`Read`], or [`BufRead`].
///
/// The loading is performed in a strict validation.
pub trait AnnotationLoader<A> {
    /// Load annotation from a file path.
    fn load_from_path<P>(&self, path: P) -> Result<A, AnnotationLoadError>
    where
        P: AsRef<Path>,
    {
        self.load_from_read(File::open(path)?)
    }

    /// Load annotation from a reader.
    fn load_from_read<R>(&self, read: R) -> Result<A, AnnotationLoadError>
    where
        R: Read,
    {
        self.load_from_buf_read(BufReader::new(read))
    }

    /// Load annotation from a buffered reader.
    fn load_from_buf_read<R>(&self, read: R) -> Result<A, AnnotationLoadError>
    where
        R: BufRead;
}
