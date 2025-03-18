use std::{io::BufRead, str::FromStr};

use ontolius::TermId;

use crate::io::AnnotationLoader;

const HPOA_COLUMN_COUNT: usize = 12;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceCode {
    /// Inferred from electronic evidence.
    IEA,

    /// Traceable author statement.
    TAS,

    /// Published clinical study.
    PCS,
}

impl FromStr for EvidenceCode {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "IEA" => Ok(EvidenceCode::IEA),
            "TAS" => Ok(EvidenceCode::TAS),
            "PCS" => Ok(EvidenceCode::PCS),
            _ => anyhow::bail!("Unknown code {s:?}")
        }
        
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Sex {
    Unknown,
    Male,
    Female,
}

impl FromStr for Sex {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Male" | "MALE" | "male" => Ok(Sex::Male),
            "Female" | "FEMALE" | "female" => Ok(Sex::Female),
            "Unknown" | "UNKNOWN" | "unknown" => Ok(Sex::Unknown),
            _ => anyhow::bail!("Unknown sex {s:?}"),
        }
    }
}

/// One of P (Phenotypic abnormality),
/// I (inheritance),
/// C (onset and clinical course).
///
/// This field is mandatory; cardinality 1
enum Aspect {
    /// Phenotypic abnormality.
    Phenotype,
    /// Inheritance.
    Inheritance,
    /// Onset and clinical course.
    ClinicalModifier,
    /// Modifier.
    Modifier,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotationReference {
    term_id: TermId,
    evidence_code: EvidenceCode,
}

struct HpoAnnotationLine<'a> {
    disease_id: &'a str,
    disease_name: &'a str,
    is_negated: bool,
    phenotype_term_id: &'a str,
    annotation_references: Vec<AnnotationReference>,
    onset: Option<&'a str>,
    frequency: &'a str,
    sex: Option<Sex>,
    modifiers: Vec<TermId>,
    aspect: Aspect,
    curators: Vec<&'a str>,
}

fn parse_line(line: &str) -> Option<HpoAnnotationLine<'_>> {
    /*
    fields = line.strip().split('\t')

    disease_id = fields[0]
    disease_name = fields[1]
    is_negated = fields[2].upper() == 'NOT'
    phenotype_id = fields[3]
    evidence_code = EvidenceCode.parse(fields[5])
    annotation_references = [AnnotationReference(TermId.from_curie(term_id), evidence_code)
                             for term_id
                             in filter(lambda t: t and not t.isspace(), fields[4].split(';'))]
    # TODO - implement parsing of temporal data
    onset = None

    frequency = fields[7]
    sex = Sex.parse(fields[8])

    modifiers = [TermId.from_curie(term_id)
                 for term_id
                 in filter(lambda t: t and not t.isspace(), fields[9].split(';'))]
    aspect = Aspect.parse(fields[10])
    curators = [curator.strip() for curator in fields[11].split(';')]

    return HpoAnnotationLine(disease_id, disease_name, is_negated,
                             phenotype_id,
                             annotation_references, onset, frequency,
                             sex, modifiers, aspect, curators)
     */
    let fields: Vec<_> = line.trim().split("\t").collect();
    if fields.len() != HPOA_COLUMN_COUNT {
        None
    } else {
        let ec = match fields[5].parse() {
            Ok(code) => code,
            Err(e) => return None,
        };
        let annotation_references = fields[4]
            .split(";")
            .filter(|&f| !f.trim().is_empty())
            .map(|curie| curie.parse().unwrap())
            .map(|tid| AnnotationReference {
                term_id: tid,
                evidence_code: ec,
            })
            .collect();
        Some(HpoAnnotationLine {
            disease_id: fields[0],
            disease_name: fields[1],
            is_negated: fields[2].eq_ignore_ascii_case("NOT"),
            phenotype_term_id: fields[3],
            annotation_references,
            onset: None, // TODO: parse
            frequency: fields[7],
            sex: fields[8].parse().ok(),
            // TODO: finalize!
            modifiers: vec![],
            aspect: Aspect::ClinicalModifier,
            curators: vec![],

        })
    }

}


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HpoAnnotation {
    term_id: TermId,
    numerator: u32,
    denominator: u32,
    references: Vec<AnnotationReference>,
    modifiers: Vec<TermId>,
}

pub struct HpoAnnotations {
    pub annotations: Vec<HpoAnnotation>,
    pub version: String,
}

pub struct HpoAnnotationLoader;

impl AnnotationLoader<HpoAnnotations> for HpoAnnotationLoader {
    fn load_from_buf_read<R>(&self, read: R) -> anyhow::Result<HpoAnnotations>
    where
        R: BufRead,
    {
        todo!()
    }
}
