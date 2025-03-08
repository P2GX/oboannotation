use flate2::bufread::GzDecoder;
use num::Integer;
use ontolius::base::TermId;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

use crate::go_relation::GoTermRelation;


#[derive(Clone, Debug, PartialEq, Serialize)]
enum EviCode {
    EXP,           // inferred from experiment
    HTP,           // Inferred from High Throughput Experiment
    PHYLO,         // Phylogenetically inferred annotations
    COMPUTATIONAL, // computational analysis evidence codes i
    AUTHOR,        // Author statement evidence
    IC,            // Curator statement
    ND,            // No biological Data available
    IEA,           // Inferred from Electronic Annotation (IEA)
}

impl FromStr for EviCode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "EXP" | "IDA" | "IPI" | "IMP" | "IGI" | "IEP" => Ok(EviCode::EXP),
            "HTP" | "HDA" | "HMP" | "HGI" | "HEP" => Ok(EviCode::HTP),
            "IBA" | "IBD" | "IKR" | "IRD" => Ok(EviCode::PHYLO),
            "ISS" | "ISO" | "ISA" | "ISM" | "IGC" | "RCA" => Ok(EviCode::COMPUTATIONAL),
            "TAS" | "NAS" => Ok(EviCode::AUTHOR),
            "IC" => Ok(EviCode::IC),
            "ND" => Ok(EviCode::ND),
            "IEA" => Ok(EviCode::IEA),
            _ => Err(format!(
                "Did not recognize '{}' as EvidenceCode.",
                s
            )),
        }
    }
}
#[derive(Serialize, Clone)]
enum Aspect {
    F,
    P,
    C,
}

impl FromStr for Aspect {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "F" => Ok(Aspect::F),
            "P" => Ok(Aspect::P),
            "C" => Ok(Aspect::C),
            _ => Err(format!(
                "Did not recognize '{}' as Aspect.",
                s
            )),
        }
    }
}



/// A Gene Ontology Annotation, corresponding to one line of the GOA file
/// 
/// We only store a subset of the information that is important for the analysis
#[derive(Clone)]
struct GoAnnot {
    gene_product_id: TermId,
    gene_product_symbol: String,
    relation: GoTermRelation,
    go_id: TermId,
    evidence_code: EviCode,
    aspect: Aspect,
}

impl GoAnnot {
    pub fn new<T: Into<String>>(
        term: TermId,
        symbol: T,
        relation: GoTermRelation,
        gene_ontology_id: TermId,
        evicode: EviCode,
        aspect: Aspect,
    ) -> Self {
        GoAnnot {
            gene_product_id: term,
            gene_product_symbol: symbol.into(),
            relation: relation,
            go_id: gene_ontology_id,
            evidence_code: evicode,
            aspect: aspect,
        }
    }
    /// A negated annotation in GO means the gene product does not have the annotation in question
    pub fn is_negated(&self) -> bool {
        return self.relation == GoTermRelation::NegatedAnnotation;
    }
}


/// To be used for serialization to display the most interesting characteristics of the annotation as a table
#[derive(Serialize)]
struct AnnotationStat {
    key: String,
    value: String,
}

impl AnnotationStat {
    pub fn from_string(item: &str, val: &str) -> Self {
        AnnotationStat {
            key: item.to_string(),
            value: val.to_string(),
        }
    }

    pub fn from_int<T: Integer + std::fmt::Display>(item: &str, val: T) -> Self {
        AnnotationStat {
            key: item.to_string(),
            value: format!("{}", val),
        }
    }
}

fn annotation_descriptive_stats(go_annots: &Vec<GoAnnot>) -> Vec<AnnotationStat> {
    let mut annots = Vec::new();
    let annot_count = go_annots.len();
    annots.push(AnnotationStat::from_int("Total annotations", annot_count));
    let unique_symbols: HashSet<_> = go_annots
        .iter()
        .map(|annot| &annot.gene_product_symbol)
        .collect();
    annots.push(AnnotationStat::from_int("genes", unique_symbols.len()));
    // Count relation types
    let mut relation_counts = HashMap::new();

    for annot in go_annots {
        *relation_counts.entry(annot.relation.clone()).or_insert(0) += 1;
    }
    for (relation, count) in &relation_counts {
        annots.push(AnnotationStat::from_int(&relation.to_string(), *count));
    }
    annots
}



pub struct GoAnnotations {
    annotation_list: Vec<GoAnnot>,
    version: String,
}


impl GoAnnotations {

    const GOA_EXPECTED_FIELDS: usize = 17;

    pub fn new(path: &str) -> Result<Self, String> {
        let file = File::open(path).expect("Could not open goa file"); // todo better error handling
        let buf_reader = BufReader::new(file); 
        let decoder = GzDecoder::new(buf_reader); 
        let reader = BufReader::new(decoder); 
        let mut annotations = vec![];
        let mut parsed_date = false; // The GOA format has multiple entries for date-generated. We only want the first
        let mut version: String = "unknown".to_string();
        for line in reader.lines() {
            match line {
                Ok(content) => {
                    if content.starts_with("!") {
                        if content.starts_with("!date-generated: ") && ! parsed_date {
                            let date_gen = &content[("!date-generated: ".len() + 2)..];
                            version = date_gen.to_string();
                            parsed_date = true;
                        }
                    } else {
                        let goann = Self::process_annotation_line(&content);
                        match goann {
                            Ok(go_annotation) => annotations.push(go_annotation),
                            Err(e) => println!("{}", e),
                        }
                    }
                },
                Err(e) => return Err(format!("Error reading file: {}", e)),
            }
        }
        print!("Parsed {} annotations", annotations.len());
        Ok(Self {
            annotation_list: annotations,
            version: version
        })
    }



    /// Process a line in go-annotation-file-gaf-format-2.2
    fn process_annotation_line(line: &str) -> Result<GoAnnot, String> {
        let tokens: Vec<&str> = line.split('\t').collect();
        if tokens.iter().count() != Self::GOA_EXPECTED_FIELDS {
            return Err(format!(
                "GOA lines expected to have {} fields, but line had {} fields: {}",
                Self::GOA_EXPECTED_FIELDS,
                tokens.iter().count(),
                line
            ));
        }
        let gene_product_id: TermId = (tokens[0], tokens[1]).into(); // Ontolius TermId
        let symbol = tokens[2];
        let relation = GoTermRelation::from_str(tokens[3])?; // return on error immediately
        let go_id = TermId::from_str(tokens[4]).unwrap(); // return on error immediately
        let evidence = EviCode::from_str(tokens[6])?; // return on error immediately
        let aspect = Aspect::from_str(tokens[8])?; // return on error immediately
        Ok(GoAnnot::new(
            gene_product_id,
            symbol,
            relation,
            go_id,
            evidence,
            aspect,
        ))
    }


    pub fn get_annotation_stats(&self) -> Vec<AnnotationStat> {
        let mut annotation_stats: Vec<AnnotationStat> = vec![];
        annotation_stats.push(AnnotationStat::from_string("version", &self.version));
        let annot_count = annotation_stats.len();
        annotation_stats.push(AnnotationStat::from_int("Total annotations", annot_count));
        let unique_symbols: HashSet<_> = self.annotation_list
            .iter()
            .map(|annot| &annot.gene_product_symbol)
            .collect();
        annotation_stats.push(AnnotationStat::from_int("genes", unique_symbols.len()));
        // Count relation types
        let mut relation_counts = HashMap::new();
    
        for annot in &self.annotation_list {
            *relation_counts.entry(annot.relation.clone()).or_insert(0) += 1;
        }
        for (relation, count) in &relation_counts {
            annotation_stats.push(AnnotationStat::from_int(&relation.to_string(), *count));
        }
        
        annotation_stats
    }

    pub fn get_annotation_statistics_json(&self) -> Result<String, String> {
        let annot_stats = self.get_annotation_stats();
        serde_json::to_string(&annot_stats).map_err(|e| format!("Serialization error: {}", e))
    }

    pub fn get_annotation_map(&self) -> HashMap<String, HashSet<TermId>> {
        let mut annot_map = HashMap::new();
        for annot in &self.annotation_list {
            if annot.is_negated() {
                continue;
            }
            let symbol = annot.gene_product_symbol.clone();
            let tid = annot.go_id.clone();
            // Insert into HashMap, creating a new HashSet if necessary
            annot_map
                .entry(symbol) 
                .or_insert_with(|| HashSet::new())  
                .insert(tid);  //
        }

        annot_map
    }
    
}




#[cfg(test)]
mod test {
    use std::assert_eq;

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
