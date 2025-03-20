use std::{io::BufRead, str::FromStr};

use anyhow::{Context, bail};
use ontolius::TermId;
use regex::Regex;
use thiserror::Error;

use crate::io::{AnnotationLoadError, AnnotationLoader, ValidationIssue};

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
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "IEA" => Ok(EvidenceCode::IEA),
            "TAS" => Ok(EvidenceCode::TAS),
            "PCS" => Ok(EvidenceCode::PCS),
            _ => Err(format!("Unknown code {s:?}")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sex {
    Unknown,
    Male,
    Female,
}

impl FromStr for Sex {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Male" | "MALE" | "male" => Ok(Sex::Male),
            "Female" | "FEMALE" | "female" => Ok(Sex::Female),
            "Unknown" | "UNKNOWN" | "unknown" => Ok(Sex::Unknown),
            _ => Err(format!("Unknown sex {s:?}")),
        }
    }
}

/// One of:
/// * `P` (Phenotypic abnormality)
/// * `I` (Mode of inheritance)
/// * `C` (Clinical course)
/// * `M` (Modifier)
/// * `H` (Past medical history)
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Aspect {
    /// Phenotypic abnormality.
    Phenotype,
    /// Inheritance.
    Inheritance,
    /// Onset and clinical course.
    ClinicalModifier,
    /// Modifier.
    Modifier,
    /// Past medical history.
    PastMedicalHistory,
}

impl TryFrom<char> for Aspect {
    type Error = String;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'P' | 'p' => Ok(Aspect::Phenotype),
            'I' | 'i' => Ok(Aspect::Inheritance),
            'C' | 'c' => Ok(Aspect::ClinicalModifier),
            'M' | 'm' => Ok(Aspect::Modifier),
            'H' | 'h' => Ok(Aspect::PastMedicalHistory),
            _ => Err(format!("Unknown aspect code: `{value}`")),
        }
    }
}

impl FromStr for Aspect {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.chars().count() {
            1 => Aspect::try_from(
                s.chars()
                    .next()
                    .expect("We just checked for presence of a single character"),
            ),
            _ => Err(format!("Unknown aspect code: `{s}`")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotationReference {
    term_id: TermId,
    evidence_code: EvidenceCode,
}

/// Frequency of a phenotypic abnormality.
#[derive(Debug, Clone, PartialEq)]
pub enum Frequency {
    /// Frequency data as [`TermId`], a member of HPO's
    /// [Frequency \[HP:0040279\]](https://hpo.jax.org/browse/term/HP:0040279) submodule.
    TermId(TermId),
    /// A count of patients affected within a cohort.
    ///
    /// For instance, `7/13` would indicate that `7 of the `13` patients with the specified disease
    /// were found to have the phenotypic abnormality referred to by the HPO term in question
    /// in the study referred to by the DB reference.
    Ratio {
        /// Count of individuals with the annotation.
        numerator: u32,
        /// The total number of the individuals investigated
        /// for the presence of the annotation.
        denominator: u32,
    },
    /// A percentage value such as 17%, again referring to the percentage of patients
    /// found to have the phenotypic abnormality referred to by the HPO term in question
    /// in the study referred to by the DB reference.
    ///
    /// If possible, the 7/13 format (see [`Frequency::Ratio`]) is preferred over the percentage format
    /// if the exact data is available.
    ///
    /// Should be a value in range of `[0..100]`. E.g. `7.` to represent 7%.
    Frequency(f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct HpoAnnotationLine {
    pub disease_id: TermId,
    pub disease_name: String,
    pub is_negated: bool,
    pub phenotype_term_id: TermId,
    pub annotation_references: Vec<AnnotationReference>,
    pub onset: Option<TermId>,
    pub frequency: Option<Frequency>,
    pub sex: Option<Sex>,
    pub modifiers: Vec<TermId>,
    pub aspect: Aspect,
    pub curators: Vec<String>,
}

/// The HPO annotation corpus with entries parsed to the highest level possible
/// without making any wild assumptions on top of the HPOA format description.
#[derive(Debug, Clone)]
pub struct HpoAnnotationLines {
    pub lines: Vec<HpoAnnotationLine>,
    pub version: String,
}

#[derive(Debug, Error)]
enum FrequencyParseError<E> {
    EmptyVal,
    Error(E),
}

struct HpoAnnotationLineParser {
    ratio_pt: Regex,
    frequency_pt: Regex,
}

impl Default for HpoAnnotationLineParser {
    fn default() -> Self {
        Self {
            ratio_pt: Regex::new(r"(?<numerator>\d+)/(?<denominator>\d+)")
                .expect("The ratio pattern should be well formatted"),
            frequency_pt: Regex::new(r"(?<frequency>\d+(\.\d*)?)%")
                .expect("The frequency pattern should be well formatted"),
        }
    }
}

impl HpoAnnotationLineParser {
    fn parse_line(&self, line: &str) -> anyhow::Result<HpoAnnotationLine> {
        let fields: Vec<_> = line.trim().split("\t").collect();
        if fields.len() != HPOA_COLUMN_COUNT {
            bail!("Cannot parse a record with {}!=13 fields", fields.len())
        } else {
            // Annotation references
            let evidence_code: EvidenceCode = match fields[5].parse() {
                Ok(ec) => ec,
                Err(e) => bail!(e),
            };

            let mut annotation_references = vec![];
            for (i, f) in fields[4].split(";").enumerate() {
                if !f.trim().is_empty() {
                    let term_id = f
                        .parse()
                        .context(format!("Parsing annotation reference #{i}"))?;
                    annotation_references.push(AnnotationReference {
                        term_id,
                        evidence_code,
                    });
                }
            }

            // Frequency
            let frequency = match self.parse_frequency(fields[7]) {
                Ok(frequency) => Some(frequency),
                Err(FrequencyParseError::EmptyVal) => None,
                Err(FrequencyParseError::Error(e)) => bail!(e),
            };

            // Modifiers
            let mut modifiers = vec![];
            for (i, f) in fields[9].split(";").enumerate() {
                if !f.trim().is_empty() {
                    modifiers.push(f.parse().context(format!("Parsing modifier #{i:?}"))?);
                }
            }

            // Aspect
            let aspect: Aspect = match fields[10].parse() {
                Ok(a) => a,
                Err(e) => bail!(e),
            };

            // The rest
            Ok(HpoAnnotationLine {
                disease_id: fields[0].parse()?,
                disease_name: fields[1].to_string(),
                is_negated: fields[2].eq_ignore_ascii_case("NOT"),
                phenotype_term_id: fields[3].parse()?,
                annotation_references,
                onset: fields[6].parse().ok(),
                frequency,
                sex: fields[8].parse().ok(),
                modifiers,
                aspect,
                curators: fields[11]
                    .split(";")
                    .map(|f| f.trim().to_string())
                    .collect(),
            })
        }
    }

    fn parse_frequency(&self, val: &str) -> Result<Frequency, FrequencyParseError<anyhow::Error>> {
        if val.is_empty() {
            Err(FrequencyParseError::EmptyVal)
        } else if let Some(cap) = self.ratio_pt.captures(val) {
            Ok(Frequency::Ratio {
                numerator: cap["numerator"]
                    .parse()
                    .expect("Regexp should ensure that numerator is parsable into a `u32` value"),
                denominator: cap["denominator"]
                    .parse()
                    .expect("Regexp should ensure that denominator is parsable into a `u32` value"),
            })
        } else if let Some(cap) = self.frequency_pt.captures(val) {
            let frequency: f64 = cap["frequency"]
                .parse()
                .expect("Regexp pattern should ensure that frequency is parsable into f64");
            if (0. ..=100.).contains(&frequency) {
                Ok(Frequency::Frequency(frequency))
            } else {
                Err(FrequencyParseError::Error(anyhow::anyhow!(
                    "Frequency not in range [0, 100]"
                )))
            }
        } else {
            // Fall back to TermId
            match val.parse().map(Frequency::TermId) {
                Ok(frequency) => Ok(frequency),
                Err(e) => Err(FrequencyParseError::Error(e)),
            }
        }
    }
}

/// Loader for HPO annotations.
///
/// ## Examples
///
/// ```
/// use oboannotation::hpo::HpoAnnotationLoader;
/// use oboannotation::io::AnnotationLoader;
///
/// let loader = HpoAnnotationLoader::default();
///
/// let data = loader.load_from_path("data/phenotype.real-shortlist.hpoa").expect("The example data should be well formatted");
///
/// // Loaded HPO annotations version `2023-04-05` ...
/// assert_eq!(data.version.as_str(), "2023-04-05");
///
/// // ... consisting of 86 lines.
/// assert_eq!(data.lines.len(), 86);
/// ```
pub struct HpoAnnotationLoader {
    line_parser: HpoAnnotationLineParser,
    version_pt: Regex,
}

impl Default for HpoAnnotationLoader {
    fn default() -> Self {
        Self {
            line_parser: Default::default(),
            version_pt: Regex::new(r"^#(date|version): (?<version>[\w-]+)\w?$")
                .expect("Default pattern should be valid"),
        }
    }
}

impl AnnotationLoader<HpoAnnotationLines> for HpoAnnotationLoader {
    fn load_from_buf_read<R>(&self, mut read: R) -> Result<HpoAnnotationLines, AnnotationLoadError>
    where
        R: BufRead,
    {
        let mut lines = vec![];
        let mut errors = vec![];
        let mut version = None;

        let mut line = String::new();
        let mut expecting_header = true;
        let mut line_number = 0usize;

        loop {
            match read.read_line(&mut line) {
                Ok(n) => {
                    if n == 0 {
                        break; // EOF was reached.
                    }

                    if expecting_header {
                        // Header
                        if line.starts_with("#DatabaseID") || line.starts_with("database_id") {
                            expecting_header = false;
                        } else if let Some(caps) = self.version_pt.captures(line.trim()) {
                            version = Some(caps["version"].to_string());
                        }
                    } else {
                        // Data
                        match self.line_parser.parse_line(&line) {
                            Ok(hal) => lines.push(hal),
                            Err(e) => errors.push(ValidationIssue::new(line_number, e.to_string())),
                        }
                    }
                }
                Err(e) => return Err(e.into()),
            };
            line.clear();
            line_number += 1;
        }

        if !errors.is_empty() {
            Err(AnnotationLoadError::ValidationError(errors))
        } else if let Some(version) = version {
            Ok(HpoAnnotationLines { lines, version })
        } else {
            Err(AnnotationLoadError::Error("Missing version".into()))
        }
    }
}
