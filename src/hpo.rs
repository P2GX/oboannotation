//! Types and I/O for working with HPO Annotations.

use std::{str::FromStr, sync::LazyLock};

use ontolius::TermId;
use regex::Regex;
use thiserror::Error;

/// Evidence codes used in HPO.
///
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum EvidenceCode {
    /// Inferred from electronic evidence.
    IEA,

    /// Traceable author statement.
    TAS,

    /// Published clinical study.
    PCS,
}

impl FromStr for EvidenceCode {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "IEA" => Ok(EvidenceCode::IEA),
            "TAS" => Ok(EvidenceCode::TAS),
            "PCS" => Ok(EvidenceCode::PCS),
            _ => Err("Unknown evidence code"),
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
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Male" | "MALE" | "male" => Ok(Sex::Male),
            "Female" | "FEMALE" | "female" => Ok(Sex::Female),
            "Unknown" | "UNKNOWN" | "unknown" => Ok(Sex::Unknown),
            _ => Err("Unknown sex"),
        }
    }
}

/// Aspect corresponds to one of the following values:
///
/// We reserve the right to add more enum variants.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
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

/// Aspect can be parsed from a `char` (case insensitive):
///
/// * Phenotypic abnormality: `P`
/// * Mode of inheritance: `I`
/// * Clinical course: `C`
/// * Modifier: `M`
/// * Past medical history: `H`
///
impl TryFrom<char> for Aspect {
    type Error = &'static str;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'P' | 'p' => Ok(Aspect::Phenotype),
            'I' | 'i' => Ok(Aspect::Inheritance),
            'C' | 'c' => Ok(Aspect::ClinicalModifier),
            'M' | 'm' => Ok(Aspect::Modifier),
            'H' | 'h' => Ok(Aspect::PastMedicalHistory),
            _ => Err("Unknown aspect code"),
        }
    }
}

/// Aspect can be parsed from a `&str` (case insensitive).
///
/// See [`Aspect::try_from<char>`] for more details.
///
impl FromStr for Aspect {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.chars().count() {
            1 => Aspect::try_from(
                s.chars()
                    .next()
                    .expect("We just checked for presence of a single character"),
            ),
            _ => Err("Unknown aspect code"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotationReference {
    term_id: TermId,
    evidence_code: EvidenceCode,
}

impl AnnotationReference {
    pub fn new(term_id: TermId, evidence_code: EvidenceCode) -> Self {
        Self {
            term_id,
            evidence_code,
        }
    }
}

/// Frequency of a phenotypic abnormality.
#[derive(Debug, Clone, PartialEq)]
pub enum Frequency {
    /// Frequency data as [`TermId`], a member of HPO's
    /// [Frequency \[HP:0040279\]](https://hpo.jax.org/browse/term/HP:0040279) submodule.
    TermId(TermId),
    /// A count of patients affected within a cohort.
    ///
    /// For instance, `7/13` would indicate that `7` of the `13` patients with the specified disease
    /// were found to have the phenotypic abnormality referred to by the HPO term in question
    /// in the study referred to by the DB reference.
    ///
    /// Note, that `1/2` and `2/4` do not represent the same information.
    /// The cohort size is `2` in the former while `4` individuals were investigated in the latter.
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

/// The possible reasons for failing to parse a frequency from a `&str`.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum FrequencyParseError {
    #[error("Empty value")]
    EmptyVal,
    #[error("Frequency not in range [0, 100]")]
    FrequencyOutOfBounds,
    #[error("Unparsable value")]
    UnparsableValue,
}

static RATIO_PT: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?<numerator>\d+)/(?<denominator>\d+)")
        .expect("The ratio pattern should be well formatted")
});

static FREQUENCY_PT: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?<frequency>-?\d+(\.\d*)?)%")
        .expect("The frequency pattern should be well formatted")
});

/// Parse a frequency string.
///
/// The parsing fails if the payload does not correspond to one of the supported input formats:
/// * ratio - e.g. `7/13` to represent 7 out of 13.
/// * percentage - e.g. 53.85%. Percent sign `%` is obligatory.
/// * term ID - e.g. `HP:0040284` for [Very rare (HP:0040284)](https://hpo.jax.org/browse/term/HP:0040284).
///
/// # Note
///
/// It is up to the user to ensure the term ID has a meaning in context of a frequency.
/// For instance, the user should check if the term
/// is a descendant of HPO's [Frequency (HP:0040279)](https://hpo.jax.org/browse/term/HP:0040279).
///
/// # Examples
///
/// ## Ratio
///
/// Parse a ratio such as `7/13`:
///
/// ```
/// use oboannotation::hpo::Frequency;
///
/// let frequency: Result<Frequency, _> = "7/13".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
/// assert!(matches!(frequency, Frequency::Ratio{numerator: 7, denominator: 13}));
/// ```
///
/// ## Percentage
///
/// Parse a percentage value, such as `53.85%`.
///
/// ```
/// use oboannotation::hpo::Frequency;
///
/// let frequency: Result<Frequency, _> = "53.85%".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
/// assert!(matches!(frequency, Frequency::Frequency(53.85)));
/// ```
///
/// ## Term id
///
/// Parse a CURIE, such as `HP:0040284`:
///
/// ```
/// use ontolius::TermId;
/// use oboannotation::hpo::Frequency;
///
/// let frequency: Result<Frequency, _> = "HP:0040284".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
/// let very_rare: TermId = "HP:0040284".parse().unwrap();
/// assert!(matches!(frequency, Frequency::TermId(very_rare)));
/// ```
///
impl FromStr for Frequency {
    type Err = FrequencyParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(FrequencyParseError::EmptyVal)
        } else if let Some(cap) = RATIO_PT.captures(s) {
            Ok(Frequency::Ratio {
                numerator: cap["numerator"]
                    .parse()
                    .expect("Regexp should ensure that numerator is parsable into a `u32` value"),
                denominator: cap["denominator"]
                    .parse()
                    .expect("Regexp should ensure that denominator is parsable into a `u32` value"),
            })
        } else if let Some(cap) = FREQUENCY_PT.captures(s) {
            let frequency: f64 = cap["frequency"]
                .parse()
                .expect("Regexp pattern should ensure that frequency is parsable into f64");
            if (0. ..=100.).contains(&frequency) {
                Ok(Frequency::Frequency(frequency))
            } else {
                Err(FrequencyParseError::FrequencyOutOfBounds)
            }
        } else {
            // Fall back to TermId
            match s.parse().map(Frequency::TermId) {
                Ok(frequency) => Ok(frequency),
                Err(_) => Err(FrequencyParseError::UnparsableValue),
            }
        }
    }
}

#[cfg(test)]
mod test_frequency {
    use crate::hpo::{Frequency, FrequencyParseError};

    #[test]
    fn from_str_empty() {
        let f: Result<Frequency, _> = "".parse();
        assert!(f.is_err());

        assert_eq!(f.unwrap_err(), FrequencyParseError::EmptyVal);
    }

    #[test]
    fn from_str_ratio() {
        let f: Result<Frequency, _> = "1/44".parse();
        assert!(f.is_ok());

        assert_eq!(
            f.unwrap(),
            Frequency::Ratio {
                numerator: 1,
                denominator: 44
            }
        );
    }

    #[test]
    fn from_str_frequency() {
        let f: Result<Frequency, _> = "54.11%".parse();
        assert!(f.is_ok());

        assert_eq!(f.unwrap(), Frequency::Frequency(54.11));
    }

    #[test]
    fn from_str_frequency_out_of_bounds() {
        assert_eq!(
            "-1.11%".parse::<Frequency>().unwrap_err(),
            FrequencyParseError::FrequencyOutOfBounds
        );
        assert_eq!(
            "100.01%".parse::<Frequency>().unwrap_err(),
            FrequencyParseError::FrequencyOutOfBounds
        );
    }

    #[test]
    fn from_str_term_id() {
        let f: Result<Frequency, _> = "HP:0040284".parse();
        assert!(f.is_ok());

        assert_eq!(f.unwrap(), Frequency::TermId("HP:0040284".parse().unwrap()));
    }

    #[test]
    fn from_str_unparsable() {
        let f: Result<Frequency, _> = "ꎯ".parse();
        assert!(f.is_err());

        assert_eq!(f.unwrap_err(), FrequencyParseError::UnparsableValue);
    }
}

/// Annotation of a disease with HPO term, including the annotation modifiers.
#[derive(Debug, Clone, PartialEq)]
pub struct HpoAnnotation {
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

/// Parse disease-phenotype annotations from HPO annotation file.
pub mod io {
    use std::io::BufRead;

    use ontolius::TermIdParseError;
    use regex::Regex;

    use crate::{
        hpo::{AnnotationReference, Aspect, EvidenceCode, FrequencyParseError, HpoAnnotation},
        io::{AnnotationLoadError, AnnotationLoader, ValidationIssue},
    };

    const HPOA_COLUMN_COUNT: usize = 12;

    /// HPOA is a tab-delimited table ...
    const DELIMITER: &str = "\t";
    /// ... with 12 columns.
    const DISEASE_ID_COL_IDX: usize = 0;
    const DISEASE_NAME_COL_IDX: usize = 1;
    const NEGATED_COL_IDX: usize = 2;
    const PHENOTYPE_ID_COL_IDX: usize = 3;
    const ANNOTATION_REFERENCES_COL_IDX: usize = 4;
    const EVIDENCE_COL_IDX: usize = 5;
    const ONSET_COL_IDX: usize = 6;
    const FREQUENCY_COL_IDX: usize = 7;
    const SEX_COL_IDX: usize = 8;
    const MODIFIERS_COL_IDX: usize = 9;
    const ASPECT_COL_IDX: usize = 10;
    const CURATORS_COL_IDX: usize = 11;

    /// The reasons for failure to parse annotation reference from a HPO annotation record.
    #[derive(Debug, thiserror::Error, PartialEq)]
    pub enum AnnotationReferenceParseError {
        #[error("Invalid annotation reference term id: #{0}")]
        InvalidTermId(usize, TermIdParseError),
    }

    #[derive(Debug)]
    pub struct HpoaError<'a> {
        fields: &'a [&'a str],
        reason: HpoaErrorReason,
    }

    impl<'a> From<(&'a [&'a str], HpoaErrorReason)> for HpoaError<'a> {
        fn from(value: (&'a [&'a str], HpoaErrorReason)) -> Self {
            Self {
                fields: value.0,
                reason: value.1,
            }
        }
    }

    impl std::fmt::Display for HpoaError<'_> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            match &self.reason {
                HpoaErrorReason::InvalidDiseaseId(e) => {
                    write!(f, "{e}: {}", self.fields[DISEASE_ID_COL_IDX])
                }
                HpoaErrorReason::InvalidPhenotypeId(_) => write!(
                    f,
                    "Invalid phenotype ID: {}",
                    self.fields[PHENOTYPE_ID_COL_IDX]
                ),
                HpoaErrorReason::InvalidEvidenceCode => {
                    write!(
                        f,
                        "Invalid evidence code: {}",
                        self.fields[EVIDENCE_COL_IDX]
                    )
                }
                HpoaErrorReason::InvalidAnnotationReference(e) => match e {
                    AnnotationReferenceParseError::InvalidTermId(i, _e) => {
                        write!(
                            f,
                            "Invalid term ID in field #{} of {}",
                            i, self.fields[ANNOTATION_REFERENCES_COL_IDX]
                        )
                    }
                },
                HpoaErrorReason::InvalidFrequency(e) => match e {
                    FrequencyParseError::EmptyVal => {
                        write!(f, "Empty value: {}", self.fields[FREQUENCY_COL_IDX])
                    }
                    FrequencyParseError::FrequencyOutOfBounds => write!(
                        f,
                        "Frequency out of bounds: {}",
                        self.fields[FREQUENCY_COL_IDX]
                    ),
                    FrequencyParseError::UnparsableValue => write!(
                        f,
                        "Unparsable frequency value: {}",
                        self.fields[FREQUENCY_COL_IDX]
                    ),
                },
                HpoaErrorReason::InvalidModifier(i, _e) => write!(
                    f,
                    "Invalid modifier ID in field #{} of {}",
                    i, self.fields[MODIFIERS_COL_IDX]
                ),
                HpoaErrorReason::InvalidAspect => {
                    write!(f, "Invalid aspect code: {}", self.fields[ASPECT_COL_IDX])
                }
                HpoaErrorReason::InvalidFieldCount(count) => {
                    write!(f, "Invalid field count {count}!={HPOA_COLUMN_COUNT}")
                }
            }
        }
    }

    impl std::error::Error for HpoaError<'_> {
        fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
            Some(&self.reason)
        }
    }

    /// The reasons for failure in parsing of HPO annotation line in [`HpoAnnotation::parse_hpoa_line`].
    #[derive(Debug, thiserror::Error, PartialEq)]
    pub enum HpoaErrorReason {
        #[error("Cannot parse a record with {0}!={HPOA_COLUMN_COUNT}")]
        InvalidFieldCount(usize),
        #[error("Invalid disease identifier: {0}")]
        InvalidDiseaseId(TermIdParseError),
        #[error("Invalid phenotype identifier: {0}")]
        InvalidPhenotypeId(TermIdParseError),
        #[error("Unknown evidence code")]
        InvalidEvidenceCode,
        #[error("Invalid annotation reference")]
        InvalidAnnotationReference(#[from] AnnotationReferenceParseError),
        #[error("Frequency parse error")]
        InvalidFrequency(#[from] FrequencyParseError),
        #[error("Modifier parse error")]
        InvalidModifier(usize, TermIdParseError),
        #[error("Unknown aspect code")]
        InvalidAspect,
    }

    /// Loader for HPO annotations.
    ///
    /// ## Examples
    ///
    /// ```
    /// use oboannotation::hpo::io::HpoAnnotationLoader;
    /// use oboannotation::io::AnnotationLoader;
    ///
    /// let loader = HpoAnnotationLoader::default();
    ///
    /// let data = loader.load_from_path("data/phenotype.real-shortlist.hpoa")
    ///              .expect("The example data should be well formatted");
    ///
    /// // Loaded HPO annotations version `2023-04-05` ...
    /// assert_eq!(data.version.as_str(), "2023-04-05");
    ///
    /// // ... consisting of 86 lines.
    /// assert_eq!(data.lines.len(), 86);
    /// ```
    pub struct HpoAnnotationLoader {
        version_pt: Regex,
    }

    impl Default for HpoAnnotationLoader {
        fn default() -> Self {
            Self {
                version_pt: Regex::new(r"^#(date|version): (?<version>[\w-]+)\w?$")
                    .expect("Default pattern should be valid"),
            }
        }
    }

    impl AnnotationLoader<HpoAnnotationLines> for HpoAnnotationLoader {
        fn load_from_buf_read<R>(
            &self,
            mut read: R,
        ) -> Result<HpoAnnotationLines, AnnotationLoadError>
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
                            match HpoAnnotation::parse_hpoa_line(&line) {
                                Ok(hal) => lines.push(hal),
                                Err(e) => {
                                    errors.push(ValidationIssue::new(line_number, &e));
                                }
                            }
                        }
                    }
                    Err(e) => return Err(e.into()),
                }
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

        fn load_from_path<P>(&self, path: P) -> Result<HpoAnnotationLines, AnnotationLoadError>
        where
            P: AsRef<std::path::Path>,
        {
            self.load_from_read(std::fs::File::open(path)?)
        }

        fn load_from_read<R>(&self, read: R) -> Result<HpoAnnotationLines, AnnotationLoadError>
        where
            R: std::io::Read,
        {
            self.load_from_buf_read(std::io::BufReader::new(read))
        }
    }

    impl HpoAnnotation {
        /// Parse HPO annotation line into `HpoAnnotation`.
        ///
        /// # Errors
        ///
        /// Parsing fails on malformed HPO annotation line.
        /// See [`HpoaErrorReason`] for the possible causes.
        ///
        pub fn parse_hpoa_line(s: &str) -> Result<Self, HpoaErrorReason> {
            let fields: Vec<_> = s.trim().split(DELIMITER).collect();
            if fields.len() == HPOA_COLUMN_COUNT {
                // Disease ID
                let disease_id = fields[DISEASE_ID_COL_IDX]
                    .parse()
                    .map_err(HpoaErrorReason::InvalidDiseaseId)?;

                // Phenotype ID
                let phenotype_term_id = fields[PHENOTYPE_ID_COL_IDX]
                    .parse()
                    .map_err(HpoaErrorReason::InvalidPhenotypeId)?;

                // Annotation references
                let mut annotation_references = vec![];
                let evidence_code: EvidenceCode = fields[EVIDENCE_COL_IDX]
                    .parse()
                    .map_err(|_e| HpoaErrorReason::InvalidEvidenceCode)?;

                for (i, f) in fields[ANNOTATION_REFERENCES_COL_IDX].split(';').enumerate() {
                    let term_id = f.parse().map_err(|e| {
                        HpoaErrorReason::InvalidAnnotationReference(
                            AnnotationReferenceParseError::InvalidTermId(i, e),
                        )
                    })?;
                    annotation_references.push(AnnotationReference {
                        term_id,
                        evidence_code: Clone::clone(&evidence_code),
                    });
                }

                // Frequency
                let frequency = match fields[FREQUENCY_COL_IDX].parse() {
                    Ok(frequency) => Some(frequency),
                    Err(e) => match e {
                        FrequencyParseError::EmptyVal => None,
                        FrequencyParseError::FrequencyOutOfBounds
                        | FrequencyParseError::UnparsableValue => {
                            return Err(HpoaErrorReason::InvalidFrequency(e));
                        }
                    },
                };

                // Modifiers
                let mut modifiers = vec![];
                for (i, f) in fields[MODIFIERS_COL_IDX].split(';').enumerate() {
                    if !f.trim().is_empty() {
                        match f.parse() {
                            Ok(t) => modifiers.push(t),
                            Err(e) => return Err(HpoaErrorReason::InvalidModifier(i, e)),
                        }
                    }
                }

                // Aspect
                let aspect: Aspect = fields[ASPECT_COL_IDX]
                    .parse()
                    .map_err(|_e| HpoaErrorReason::InvalidAspect)?;

                // The rest
                Ok(HpoAnnotation {
                    disease_id,
                    disease_name: fields[DISEASE_NAME_COL_IDX].to_string(),
                    is_negated: fields[NEGATED_COL_IDX].eq_ignore_ascii_case("NOT"),
                    phenotype_term_id,
                    annotation_references,
                    onset: fields[ONSET_COL_IDX].parse().ok(),
                    frequency,
                    sex: fields[SEX_COL_IDX].parse().ok(),
                    modifiers,
                    aspect,
                    curators: fields[CURATORS_COL_IDX]
                        .split(';')
                        .map(|f| f.trim().to_string())
                        .collect(),
                })
            } else {
                Err(HpoaErrorReason::InvalidFieldCount(fields.len()))
            }
        }
    }

    /// The HPO annotation corpus with entries parsed to the highest level possible
    /// while making no wild assumptions about the parsed data.
    #[derive(Debug, Clone)]
    pub struct HpoAnnotationLines {
        /// The HPO annotation records.
        pub lines: Vec<HpoAnnotation>,
        /// The HPOA version (e.g. `2023-04-05`)
        pub version: String,
    }

    #[cfg(test)]
    mod test_hpo_io {
        use ontolius::TermIdParseError;

        use crate::hpo::{AnnotationReference, Aspect, EvidenceCode, Frequency};

        use super::{HpoAnnotation, HpoaErrorReason};

        #[test]
        fn parse_hpoa_line_ok() {
            let line = "OMIM:154700\tMarfan syndrome\t\tHP:0001377\tPMID:28050285;PMID:33436942\tPCS\t\t29/199\t\t\tP\tHPO:probinson[2021-05-27];HPO:probinson[2021-04-01]";

            let hpo_line: Result<HpoAnnotation, _> = HpoAnnotation::parse_hpoa_line(line);

            assert!(hpo_line.is_ok());

            let hpo_line = hpo_line.unwrap();
            assert_eq!(&hpo_line.disease_id.to_string(), &"OMIM:154700");
            assert_eq!(hpo_line.disease_name.as_str(), "Marfan syndrome");
            assert_eq!(hpo_line.is_negated, false);
            assert_eq!(&hpo_line.phenotype_term_id.to_string(), &"HP:0001377");
            assert_eq!(
                &hpo_line.annotation_references,
                &[
                    AnnotationReference::new("PMID:28050285".parse().unwrap(), EvidenceCode::PCS),
                    AnnotationReference::new("PMID:33436942".parse().unwrap(), EvidenceCode::PCS)
                ]
            );
            assert!(hpo_line.onset.is_none());
            assert_eq!(hpo_line.frequency.is_some(), true);
            assert_eq!(
                &hpo_line.frequency.unwrap(),
                &Frequency::Ratio {
                    numerator: 29,
                    denominator: 199
                }
            );

            assert!(hpo_line.sex.is_none());
            assert!(hpo_line.modifiers.is_empty());
            assert_eq!(hpo_line.aspect, Aspect::Phenotype);
            assert_eq!(
                &hpo_line.curators,
                &[
                    "HPO:probinson[2021-05-27]".to_string(),
                    "HPO:probinson[2021-04-01]".to_string()
                ]
            );
        }

        #[test]
        fn parse_hpoa_line_bad_disease_id() {
            let line = "OMIM-154700\tMarfan syndrome\t\tHP:0001377\tPMID:28050285;PMID:33436942\tPCS\t\t29/199\t\t\tP\tHPO:probinson[2021-05-27];HPO:probinson[2021-04-01]";
            let hpo_line: Result<_, _> = HpoAnnotation::parse_hpoa_line(line);

            assert!(hpo_line.is_err());
            assert_eq!(
                hpo_line.unwrap_err(),
                HpoaErrorReason::InvalidDiseaseId(TermIdParseError::MissingDelimiter)
            );
        }
    }
}
