//! This library simplifies accessing annotations of biomedical ontologies.
//!
//! # Supported ontologies
//!
//! The library supports the following ontologies
//!
//! * [Human Phenotype Ontology (HPO)](https://hpo.jax.org/)
//! * [Gene Ontology (GO)](https://www.geneontology.org/)
//!
//! # Examples
//!
//! ## Load HPO annotations
//!
//! Load a toy HPO annotation file:
//!
//! ```rust
//! use oboannotation::io::AnnotationLoader;
//! use oboannotation::hpo::HpoAnnotation;
//! use oboannotation::hpo::io::{HpoAnnotationLines, HpoAnnotationLoader};
//!
//! // 👇 Replace with path to a real file 👇
//! let fpath_hpoa = "data/phenotype.real-shortlist.hpoa";
//!
//! let loader = HpoAnnotationLoader::default();
//! let data: HpoAnnotationLines = loader.load_from_path(fpath_hpoa)
//!                                  .expect("Toy HPOA should be well formatted");
//!
//! assert_eq!(data.lines.len(), 86); // Toy HPOA includes 86 records
//! ```
//!
//! See [`HpoAnnotationLines`][`crate::hpo::io::HpoAnnotationLines`] and [`HpoAnnotation`][`crate::hpo::HpoAnnotation`]
//! to learn more about the data format.
//!
pub mod go;
pub mod hpo;
pub mod io;
