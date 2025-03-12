//! Parse Gene Ontology annotations.
//!
//! Use [`GoGafAnnotationLoader`] to parse GAF file into [`GoAnnotations`].
use anyhow::bail;
use ontolius::TermId;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt::Display;
use std::io::BufRead;
use std::str::FromStr;

use crate::io::AnnotationLoader;

/// The number of columns in GO GAF file.
const GOA_EXPECTED_FIELDS: usize = 17;

#[derive(Debug)]
pub enum InputError {
    NegatedAnnotation, // we skip negated annotations
    MalformedLine(String),
    ParsingError(String), // Another error type
                          // Add other error kinds as needed
}

impl Display for InputError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            InputError::NegatedAnnotation => write!(f, "Negated annotation detected"),
            InputError::MalformedLine(ref s) => write!(f, "Malformed line: {}", s),
            InputError::ParsingError(ref s) => write!(f, "Parsing error: {}", s),
        }
    }
}

impl Error for InputError {}

/// Gene product to GO term relations
/// enables links a gene product to a Molecular Function it executes.
/// contributes to links a gene product to a Molecular Function executed by a macromolecular complex, in which the Molecular Function cannot be ascribed to an individual subunit of that complex. Only the complex subunits required for the Molecular Function are annotated to the Molecular Function term with ‘contributes to’.
/// Relations between a gene product and a Biological Process:
/// involved in links a gene product and a Biological Process in which the gene product’s Molecular Function plays an integral role.
/// acts upstream of or within links a gene product and a Biological Process when the mechanism relating the gene product’s activity to the Biological Process is not known.
/// Relations between a gene product and a Cellular Component:
/// is active in links a gene product to the cellular location in which it enables its Molecular Function.
/// located in links a gene product and the Cellular Component, specifically a cellular anatomical anatomy or virion component, in which a gene product has been detected.
/// part of links a gene product and a protein-containing complex.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum GoTermRelation {
    Enables,
    ContributesTo,
    InvolvedIn,
    ActsUpstreamOf,
    ActsWithin,
    ActsUpstreamOfOrWithin,
    ActsUpstreamOfNegativeEffect,
    ActsUpstreamOfPositiveEffect,
    ActsUpstreamOfOrWithinNegativeEffect,
    ActsUpstreamOfOrWithinPositiveEffect,
    IsActiveIn,
    LocatedIn,
    ColocalizesWith,
    PartOf,
}

impl Display for GoTermRelation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let relation_str = match self {
            GoTermRelation::Enables => "enables",
            GoTermRelation::ContributesTo => "contributes_to",
            GoTermRelation::InvolvedIn => "involved_in",
            GoTermRelation::ActsUpstreamOf => "acts_upstream_of",
            GoTermRelation::ActsWithin => "acts_within",
            GoTermRelation::IsActiveIn => "is_active_in",
            GoTermRelation::ActsUpstreamOfOrWithin => "acts_upstream_of_or_within",
            GoTermRelation::ActsUpstreamOfNegativeEffect => "acts_upstream_of_negative_effect",
            GoTermRelation::ActsUpstreamOfPositiveEffect => "acts_upstream_of_positive_effect",
            GoTermRelation::ActsUpstreamOfOrWithinNegativeEffect => {
                "acts_upstream_of_or_within_negative_effect"
            }
            GoTermRelation::ActsUpstreamOfOrWithinPositiveEffect => {
                "acts_upstream_of_or_within_positive_effect"
            }
            GoTermRelation::LocatedIn => "located_in",
            GoTermRelation::PartOf => "part_of",
            GoTermRelation::ColocalizesWith => "colocalizes_with",
        };
        write!(f, "{}", relation_str)
    }
}

impl FromStr for GoTermRelation {
    type Err = InputError;

    fn from_str(s: &str) -> Result<Self, InputError> {
        if s.starts_with("NOT") {
            return Err(InputError::NegatedAnnotation);
        }
        match s {
            "enables" => Ok(GoTermRelation::Enables),
            "contributes_to" => Ok(GoTermRelation::ContributesTo),
            "involved_in" => Ok(GoTermRelation::InvolvedIn),
            "located_in" => Ok(GoTermRelation::LocatedIn),
            "acts_upstream_of" => Ok(GoTermRelation::ActsUpstreamOf),
            "acts_within" => Ok(GoTermRelation::ActsWithin),
            "acts_upstream_of_or_within" => Ok(GoTermRelation::ActsUpstreamOfOrWithin),
            "acts_upstream_of_negative_effect" => Ok(GoTermRelation::ActsUpstreamOfNegativeEffect),
            "acts_upstream_of_positive_effect" => Ok(GoTermRelation::ActsUpstreamOfPositiveEffect),
            "acts_upstream_of_or_within_negative_effect" => {
                Ok(GoTermRelation::ActsUpstreamOfOrWithinNegativeEffect)
            }
            "acts_upstream_of_or_within_positive_effect" => {
                Ok(GoTermRelation::ActsUpstreamOfOrWithinPositiveEffect)
            }
            "is_active_in" => Ok(GoTermRelation::IsActiveIn),
            "part_of" => Ok(GoTermRelation::PartOf),
            "colocalizes_with" => Ok(GoTermRelation::ColocalizesWith),
            _ => Err(InputError::ParsingError(format!(
                "Did not recognize '{}' as a GOA relation.",
                s
            ))),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EviCode {
    EXP,           // Inferred from experiment
    HTP,           // Inferred from High Throughput Experiment
    PHYLO,         // Phylogenetically inferred annotations
    COMPUTATIONAL, // Computational analysis evidence codes i
    AUTHOR,        // Author statement evidence
    IC,            // Curator statement
    ND,            // No biological Data available
    IEA,           // Inferred from Electronic Annotation (IEA)
}

impl FromStr for EviCode {
    type Err = InputError;

    fn from_str(s: &str) -> Result<Self, InputError> {
        match s {
            "EXP" | "IDA" | "IPI" | "IMP" | "IGI" | "IEP" => Ok(EviCode::EXP),
            "HTP" | "HDA" | "HMP" | "HGI" | "HEP" => Ok(EviCode::HTP),
            "IBA" | "IBD" | "IKR" | "IRD" => Ok(EviCode::PHYLO),
            "ISS" | "ISO" | "ISA" | "ISM" | "IGC" | "RCA" => Ok(EviCode::COMPUTATIONAL),
            "TAS" | "NAS" => Ok(EviCode::AUTHOR),
            "IC" => Ok(EviCode::IC),
            "ND" => Ok(EviCode::ND),
            "IEA" => Ok(EviCode::IEA),
            _ => Err(InputError::ParsingError(format!(
                "Did not recognize '{}' as EvidenceCode.",
                s
            ))),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Aspect {
    /// Molecular function.
    F,
    /// Biological process.
    P,
    /// Cellular component.
    C,
}

impl FromStr for Aspect {
    type Err = InputError;

    fn from_str(s: &str) -> Result<Self, InputError> {
        match s {
            "F" => Ok(Aspect::F),
            "P" => Ok(Aspect::P),
            "C" => Ok(Aspect::C),
            _ => Err(InputError::ParsingError(format!(
                "Did not recognize '{}' as Aspect.",
                s
            ))),
        }
    }
}

/// A Gene Ontology Annotation, corresponding to one line of the GOA file
///
/// We only store a subset of the information that is important for the analysis
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GoAnnotation {
    gene_product_id: TermId,
    gene_product_symbol: String,
    relation: GoTermRelation,
    gene_ontology_id: TermId,
    evidence_code: EviCode,
    aspect: Aspect,
}

impl GoAnnotation {
    pub fn new<T>(
        term: TermId,
        symbol: T,
        relation: GoTermRelation,
        gene_ontology_id: TermId,
        evidence_code: EviCode,
        aspect: Aspect,
    ) -> Self
    where
        T: ToString,
    {
        GoAnnotation {
            gene_product_id: term,
            gene_product_symbol: symbol.to_string(),
            relation,
            gene_ontology_id,
            evidence_code,
            aspect,
        }
    }
}

/// A container with all GO annotations, as parsed from the annotation data file.
pub struct GoAnnotations {
    pub annotations: Vec<GoAnnotation>,
    pub version: String,
    pub negated_annotation_count: usize,
}

/// To be used for serialization to display the most interesting characteristics of the annotation as a table
// #[derive(Debug, Clone, PartialEq, Eq)]
// #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
// struct AnnotationStat {
//     key: String,
//     value: String,
// }

// impl AnnotationStat {
//     pub fn from_string(item: &str, val: &str) -> Self {
//         AnnotationStat {
//             key: item.to_string(),
//             value: val.to_string(),
//         }
//     }

//     pub fn from_int<T>(item: &str, val: T) -> Self
//     where
//         T: ToString,
//     {
//         AnnotationStat {
//             key: item.to_string(),
//             value: val.to_string(),
//         }
//     }
// }

// fn annotation_descriptive_stats(go_annots: &[GoAnnotation]) -> Vec<AnnotationStat> {
//     let mut annots = Vec::new();
//     let annot_count = go_annots.len();
//     annots.push(AnnotationStat::from_int("Total annotations", annot_count));
//     let unique_symbols: HashSet<_> = go_annots
//         .iter()
//         .map(|annot| &annot.gene_product_symbol)
//         .collect();
//     annots.push(AnnotationStat::from_int("genes", unique_symbols.len()));
//     // Count relation types
//     let mut relation_counts = HashMap::new();

//     for annot in go_annots {
//         *relation_counts.entry(annot.relation.clone()).or_insert(0) += 1;
//     }
//     for (relation, count) in &relation_counts {
//         annots.push(AnnotationStat::from_int(&relation.to_string(), *count));
//     }
//     annots
// }

pub struct GoGafAnnotationLoader;

impl AnnotationLoader<GoAnnotations> for GoGafAnnotationLoader {
    fn load_from_buf_read<R>(&self, read: R) -> anyhow::Result<GoAnnotations>
    where
        R: BufRead,
    {
        let mut annotations = vec![];
        let mut negated_annotation_count = 0;
        let mut parsed_date = false; // The GOA format has multiple entries for date-generated. We only want the first
        let mut version = "UNKNOWN".to_string();
        for line in read.lines() {
            match line {
                Ok(line) => {
                    if line.starts_with("!") {
                        if line.starts_with("!date-generated: ") && !parsed_date {
                            let date_gen = &line[("!date-generated: ".len() + 2)..];
                            parsed_date = true;
                            version = date_gen.to_string();
                        }
                    } else {
                        match parse_annotation_line(&line) {
                            Ok(ann) => annotations.push(ann),
                            Err(InputError::NegatedAnnotation) => negated_annotation_count += 1,
                            Err(e) => bail!(e),
                        }
                    }
                }
                Err(e) => bail!(e),
            }
        }

        Ok(GoAnnotations {
            annotations,
            negated_annotation_count,
            version,
        })
    }
}

/// Process a line in go-annotation-file-gaf-format-2.2
fn parse_annotation_line(line: &str) -> Result<GoAnnotation, InputError> {
    let tokens: Vec<&str> = line.split('\t').collect();
    if tokens.len() != GOA_EXPECTED_FIELDS {
        return Err(InputError::MalformedLine(format!(
            "GOA lines expected to have {} fields, but line had {} fields: {}",
            GOA_EXPECTED_FIELDS,
            tokens.len(),
            line
        )));
    }
    let gene_product_id = TermId::from((tokens[0], tokens[1]));
    let symbol = tokens[2];
    let relation = GoTermRelation::from_str(tokens[3])?; // return on error immediately
    let go_id = match tokens[4].parse() {
        Ok(t) => t,
        Err(e) => return Err(InputError::ParsingError(format!("Error: {e:?}"))),
    }; // return on error immediately
    let evidence = EviCode::from_str(tokens[6])?; // return on error immediately
    let aspect = Aspect::from_str(tokens[8])?; // return on error immediately
    Ok(GoAnnotation::new(
        gene_product_id,
        symbol,
        relation,
        go_id,
        evidence,
        aspect,
    ))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_valid_evidence_codes() {
        let tests = vec![
            ("EXP", EviCode::EXP),
            ("IDA", EviCode::EXP),
            ("IPI", EviCode::EXP),
            ("IMP", EviCode::EXP),
            ("IEP", EviCode::EXP),
            ("HTP", EviCode::HTP),
            ("HDA", EviCode::HTP),
            ("HMP", EviCode::HTP),
            ("HGI", EviCode::HTP),
            ("HEP", EviCode::HTP),
            ("IBA", EviCode::PHYLO),
            ("IBD", EviCode::PHYLO),
            ("IKR", EviCode::PHYLO),
            ("IRD", EviCode::PHYLO),
            ("ISS", EviCode::COMPUTATIONAL),
            ("ISO", EviCode::COMPUTATIONAL),
            ("ISA", EviCode::COMPUTATIONAL),
            ("ISM", EviCode::COMPUTATIONAL),
            ("ISS", EviCode::COMPUTATIONAL),
            ("ISS", EviCode::COMPUTATIONAL),
            ("TAS", EviCode::AUTHOR),
            ("NAS", EviCode::AUTHOR),
            ("IC", EviCode::IC),
            ("ND", EviCode::ND),
            ("IEA", EviCode::IEA),
        ];
        for test in tests {
            let ecode = EviCode::from_str(test.0);
            assert!(ecode.is_ok());
            assert_eq!(ecode.unwrap(), test.1);
        }
    }

    /// Make sure we get an error with an invalid evidence code
    #[test]
    fn test_invalid_evidence_code() {
        let ecode = EviCode::from_str("something");
        assert!(ecode.is_err());
        match ecode {
            Err(e) => {
                assert_eq!(
                    e.to_string(),
                    "Parsing error: Did not recognize 'something' as EvidenceCode.".to_string()
                );
            }
            Ok(_) => panic!("Expected an error, but got Ok."),
        }
    }
}
