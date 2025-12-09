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
//! use oboannotation::hpo::{HpoAnnotation, HpoAnnotations};
//!
//! // 👇 Replace with path to a real file 👇
//! let fpath_hpoa = "data/phenotype.real-shortlist.hpoa";
//! let data: HpoAnnotations = HpoAnnotations::load_from_path(fpath_hpoa)
//!                                  .expect("Toy HPOA should be well formatted");
//!
//! assert_eq!(data.annotations().len(), 86); // Toy HPOA includes 86 records
//! ```
//!
//! See [`HpoAnnotations`][`crate::hpo::HpoAnnotations`] and [`HpoAnnotation`][`crate::hpo::HpoAnnotation`]
//! to learn more about the data format.
//!
//! ## Load GO annotations
//!
//! Load a toy GO annotation file containing annotations of *FBN1* and *SURF1* genes for *Homo sapiens*.
//!
//! ```rust
//! use oboannotation::io::AnnotationLoader;
//! use oboannotation::go::{GoAnnotations};
//!
//! // 👇 Replace with path to a real file 👇
//! let fpath_goa = "data/goa_human.SURF1_FBN1.gaf";
//!
//! let data: GoAnnotations = GoAnnotations::load_from_path(fpath_goa)
//!                             .expect("Toy GO annotations should be well formatted");
//!
//! assert_eq!(data.annotations.len(), 156); // Toy GO annotations include 156 records
//! ```
//! See [`GoAnnotations`][`crate::go::GoAnnotations`] and [`GoAnnotation`][`crate::go::GoAnnotation`]
//! to learn more about the data format.
//!
pub mod format;
pub mod go;
pub mod hpo;
pub mod io;
