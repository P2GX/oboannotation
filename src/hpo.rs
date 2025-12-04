//! Types and I/O for working with HPO Annotations.

use ontolius::TermId;
use regex::Regex;
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::marker::PhantomData;
use std::{str::FromStr, sync::LazyLock};

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
    ClinicalCourse,
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
            'C' | 'c' => Ok(Aspect::ClinicalCourse),
            'M' | 'm' => Ok(Aspect::Modifier),
            'H' | 'h' => Ok(Aspect::PastMedicalHistory),
            _ => Err("Unknown aspect code"),
        }
    }
}

/// Aspect can be parsed from a `&str` (case-insensitive).
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

/// The ways to encode frequency of an HPO annotation.
///
/// The annotation frequency can be in one of the following forms:
/// * [`FrequencyData::TermId`] - a term ID (e.g. [Obligate \[HP:0040280\]](https://hpo.jax.org/browse/term/HP:0040280))
/// * [`FrequencyData::Ratio`] - an `n` over `m` (e.g. `7` out of `13` investigated individuals presented with the feature)
/// * [`FrequencyData::Percentage`] - a percentage (e.g. 13% of the individuals presented with the feature)
///
/// The [`FrequencyData::Ratio`] is preferred over the other forms.
#[derive(Debug, Clone, PartialEq)]
pub enum FrequencyData {
    /*
    Notes
    `FrequencyData` cannot implement `Eq` because `FrequencyData::Percentage` can be created with a `f64:NAN`.
    Nor can `Hash` be implemented due to `f64`.
     */
    /// Frequency data as [`TermId`], a member of HPO's
    /// [Frequency \[HP:0040279\]](https://hpo.jax.org/browse/term/HP:0040279) submodule.
    TermId(TermId),
    /// A count of patients affected within a cohort.
    ///
    /// For instance, `7/13` would indicate that `7` (n) of the `13` (m) patients with the specified disease
    /// were found to have the phenotypic abnormality referred to by the HPO term in question
    /// in the study referred to by the DB reference.
    ///
    /// Note, that `1/2` and `2/4` do not represent the same information.
    /// The cohort size is `2` in the former while `4` individuals were investigated in the latter.
    ///
    /// `n` SHOULD be less than or equal to `m`.
    Ratio {
        /// Count of individuals with the annotation.
        n: u32,
        /// The total number of the individuals investigated
        /// for the presence of the annotation.
        m: u32,
    },
    /// A percentage value such as 17%, again referring to the percentage of patients
    /// found to have the phenotypic abnormality referred to by the HPO term in question
    /// in the study referred to by the DB reference.
    ///
    /// If possible, the 7/13 format (see [`FrequencyData::Ratio`]) is preferred over the percentage format
    /// if the exact data is available.
    ///
    /// The percentage SHOULD be a value in range of `[0..100]` (e.g. `7.` to represent 7%.).
    Percentage(f64),
}

impl From<TermId> for FrequencyData {
    fn from(value: TermId) -> Self {
        Self::TermId(value)
    }
}

/// Frequency of a phenotypic abnormality.
///
/// The frequency is guaranteed to meet the following invariants:
///
/// * in [`FrequencyData::Ratio`], `n` is less than or equal to `m`
/// * in [`FrequencyData::Percentage`], the value is in the range of \[0, 100\]
#[derive(Debug, Clone, PartialEq)]
pub struct Frequency(FrequencyData);

impl Frequency {
    /// Get the frequency data.
    pub fn data(&self) -> &FrequencyData {
        &self.0
    }

    /// Create `Frequency` from `n` over `m`.
    ///
    ///
    /// # Example
    ///
    /// Create frequency for annotation that was observed in `4` out of `8` tested individuals:
    ///
    /// ```
    /// use oboannotation::hpo::{Frequency, FrequencyData};
    ///
    /// let freq: Result<Frequency, _> = Frequency::from_ratio(4u8, 8u8);
    ///
    /// assert!(freq.is_ok());
    /// assert_eq!(freq.unwrap().data(), &FrequencyData::Ratio {n: 4, m: 8});
    /// ```
    ///
    /// # Errors
    ///
    /// Fails if `n` is greater than `m`:
    ///
    /// ```
    /// use oboannotation::hpo::{Frequency, FrequencyData, FrequencyParseError};
    ///
    /// let freq: Result<Frequency, FrequencyParseError> = Frequency::from_ratio(9u8, 8u8);
    ///
    /// assert!(freq.is_err());
    /// assert_eq!(freq.unwrap_err(), FrequencyParseError::NGreaterThanM);
    /// ```
    pub fn from_ratio(n: impl Into<u32>, m: impl Into<u32>) -> Result<Self, FrequencyParseError> {
        let (n, m) = (n.into(), m.into());
        if n <= m {
            Ok(Self(FrequencyData::Ratio { n, m }))
        } else {
            Err(FrequencyParseError::NGreaterThanM)
        }
    }

    /// Create `Frequency` from a percentage.
    ///
    /// # Example
    ///
    /// Create frequency for annotation that was observed in 49.5% of individuals:
    ///
    /// ```
    /// use oboannotation::hpo::{Frequency, FrequencyData};
    ///
    /// let freq: Result<Frequency, _> = Frequency::from_percentage(49.5);
    ///
    /// assert!(freq.is_ok());
    /// assert_eq!(freq.unwrap().data(), &FrequencyData::Percentage(49.5));
    /// ```
    ///
    /// # Errors
    ///
    /// Fails if the percentage is not in range of \[0, 100\]:
    ///
    /// ```
    /// use oboannotation::hpo::{Frequency, FrequencyData, FrequencyParseError};
    ///
    /// let freq: Result<Frequency, FrequencyParseError> = Frequency::from_percentage(100.01);
    ///
    /// assert!(freq.is_err());
    /// assert_eq!(freq.unwrap_err(), FrequencyParseError::PercentageOutOfBounds);
    /// ```
    pub fn from_percentage(percentage: impl Into<f64>) -> Result<Self, FrequencyParseError> {
        let percentage = percentage.into();

        if f64::is_sign_positive(percentage) && percentage <= 100.0 {
            Ok(Self(FrequencyData::Percentage(percentage)))
        } else {
            Err(FrequencyParseError::PercentageOutOfBounds)
        }
    }
}

/// Format the frequency.
///
/// # Examples
///
/// ```
/// use ontolius::TermId;
/// use oboannotation::hpo::Frequency;
///
/// let freq: Frequency = Frequency::from_ratio(7u8, 13u8).unwrap();
///
/// assert_eq!(freq.to_string().as_str(), "7/13");
/// ```
impl Display for Frequency {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            FrequencyData::TermId(t) => write!(f, "{}", t),
            FrequencyData::Ratio { n, m } => write!(f, "{n}/{m}"),
            FrequencyData::Percentage(percentage) => write!(f, "{:.1}%", percentage),
        }
    }
}

/// Convert [`FrequencyData`] to `Frequency` and validate its invariants.
impl TryFrom<FrequencyData> for Frequency {
    type Error = FrequencyParseError;

    fn try_from(value: FrequencyData) -> Result<Self, Self::Error> {
        match value {
            FrequencyData::TermId(t) => Ok(Self(FrequencyData::TermId(t))),
            FrequencyData::Ratio { n, m } => Self::from_ratio(n, m),
            FrequencyData::Percentage(val) => Self::from_percentage(val),
        }
    }
}

impl From<TermId> for Frequency {
    fn from(value: TermId) -> Self {
        Self(value.into())
    }
}

/// The possible reasons for failing to parse a frequency from a `&str`.
#[derive(Clone, Debug, thiserror::Error, PartialEq)]
pub enum FrequencyParseError {
    #[error("Empty value")]
    EmptyVal,
    #[error("`m` is greater than `n`")]
    NGreaterThanM,
    #[error("Percentage not in range [0, 100]")]
    PercentageOutOfBounds,
    #[error("Unparsable value")]
    UnparsableValue,
}

static RATIO_PT: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(?<numerator>\d+)/(?<denominator>\d+)$")
        .expect("The ratio pattern should be well formatted")
});

static FREQUENCY_PT: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(?<frequency>-?\d+(\.\d*)?)%$")
        .expect("The frequency pattern should be well formatted")
});

/// Parse a frequency string.
///
/// The parsing fails if the payload does not correspond to one of the supported input formats:
/// * ratio - `n`/`m` (e.g. `7/13`) to represent presence of a feature in `n` individuals
///   out of `m` tested for feature's presence.
/// * percentage - e.g. 53.85%. Percent sign `%` is obligatory. The value must be between \[0, 100\].
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
/// use oboannotation::hpo::{Frequency, FrequencyData};
///
/// let frequency: Result<Frequency, _> = "7/13".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
///
/// assert_eq!(frequency.data(), &FrequencyData::Ratio {n: 7, m: 13});
/// ```
///
/// Fails if `n > m`:
/// ```
/// use oboannotation::hpo::{Frequency, FrequencyParseError};
///
/// let frequency: FrequencyParseError = "13/7".parse::<Frequency>().unwrap_err();
///
/// assert_eq!(frequency, FrequencyParseError::NGreaterThanM);
/// ```
///
/// ## Percentage
///
/// Parse a percentage value, such as `53.85%`.
///
/// ```
/// use oboannotation::hpo::{Frequency, FrequencyData};
///
/// let frequency: Result<Frequency, _> = "53.85%".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
/// assert_eq!(frequency.data(), &FrequencyData::Percentage(53.85));
/// ```
///
/// ## Term id
///
/// Parse a CURIE, such as `HP:0040284`:
///
/// ```
/// use ontolius::TermId;
/// use oboannotation::hpo::{Frequency, FrequencyData};
///
/// let frequency: Result<Frequency, _> = "HP:0040284".parse();
/// assert!(frequency.is_ok());
///
/// let frequency = frequency.unwrap();
/// let very_rare: TermId = "HP:0040284".parse().unwrap();
/// assert_eq!(frequency.data(), &FrequencyData::TermId(very_rare));
/// ```
///
impl FromStr for Frequency {
    type Err = FrequencyParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(FrequencyParseError::EmptyVal)
        } else if let Some(cap) = RATIO_PT.captures(s) {
            let n: u32 = cap["numerator"]
                .parse()
                .expect("Regexp should ensure that numerator is parsable into a `u32` value");
            let m: u32 = cap["denominator"]
                .parse()
                .expect("Regexp should ensure that denominator is parsable into a `u32` value");
            Frequency::from_ratio(n, m)
        } else if let Some(cap) = FREQUENCY_PT.captures(s) {
            let frequency: f64 = cap["frequency"]
                .parse()
                .expect("Regexp pattern should ensure that frequency is parsable into f64");
            Frequency::from_percentage(frequency)
        } else {
            // Fall back to TermId
            match s.parse().map(|t| Self(FrequencyData::TermId(t))) {
                Ok(frequency) => Ok(frequency),
                Err(_) => Err(FrequencyParseError::UnparsableValue),
            }
        }
    }
}

#[cfg(test)]
mod test_frequency {
    use super::{Frequency, FrequencyData, FrequencyParseError};

    #[test]
    fn from_str_empty() {
        let f: Result<Frequency, _> = "".parse();
        assert!(f.is_err());

        assert_eq!(f.unwrap_err(), FrequencyParseError::EmptyVal);
    }

    #[test]
    fn from_str_freq_bad_ratio() {
        let f: Result<Frequency, _> = "13/7".parse();
        assert!(f.is_err());

        assert_eq!(f.unwrap_err(), FrequencyParseError::NGreaterThanM);
    }

    #[test]
    fn from_str_percentage_out_of_bounds() {
        assert_eq!(
            "-0%".parse::<Frequency>().unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
        assert_eq!(
            "100.01%".parse::<Frequency>().unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
    }

    #[test]
    fn from_str_unparsable() {
        let f: Result<Frequency, _> = "ꎯ".parse();
        assert!(f.is_err());

        assert_eq!(f.unwrap_err(), FrequencyParseError::UnparsableValue);
    }

    #[test]
    fn test_creating_from_negative_zero_percent() {
        let result = Frequency::from_percentage(-0.);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
    }

    #[test]
    fn test_cannot_create_from_nan_or_infinity() {
        assert_eq!(
            Frequency::from_percentage(f64::NAN).unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
        assert_eq!(
            Frequency::from_percentage(f64::INFINITY).unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
        assert_eq!(
            Frequency::from_percentage(f64::NEG_INFINITY).unwrap_err(),
            FrequencyParseError::PercentageOutOfBounds
        );
    }

    #[test]
    fn test_can_create_percentage_from_min_positive() {
        assert!(Frequency::from_percentage(f64::MIN_POSITIVE).is_ok());
    }

    #[test]
    fn test_we_pay_nothing_for_the_frequency_wrapper() {
        assert_eq!(size_of::<Frequency>(), size_of::<FrequencyData>());
    }
}

/// Annotation of a disease with HPO term, including the annotation modifiers, onset, sex,
/// clinical modifiers and the curators.
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

impl HpoAnnotation {
    /// Create a builder for building the [`HpoAnnotation`].
    pub fn builder() -> HpoAnnotationBuilder<Unset, Unset, Unset> {
        HpoAnnotationBuilder {
            disease: None,
            is_negated: false,
            phenotype_term_id: None,
            annotation_references: vec![],
            onset: None,
            frequency: None,
            sex: None,
            modifiers: vec![],
            aspect: None,
            curators: vec![],
            state: PhantomData,
        }
    }
}

/// A marker struct to indicate that a required field of [`HpoAnnotation] was set.
pub struct Set;

/// A marker struct to indicate that a required field of [`HpoAnnotation] was not set.
pub struct Unset;

/// A builder for [`HpoAnnotation`].
///
/// Three required fields must be set on the builder ...
/// * [`HpoAnnotationBuilder::disease`]
/// * [`HpoAnnotationBuilder::phenotype`]
/// * [`HpoAnnotationBuilder::aspect`]
///
/// ... before the annotation can be built with [`HpoAnnotationBuilder::build`].
pub struct HpoAnnotationBuilder<D, P, A> {
    // The builder has a generic parameter for each required field/condition.
    disease: Option<(TermId, String)>,
    is_negated: bool,
    phenotype_term_id: Option<TermId>,
    annotation_references: Vec<AnnotationReference>,
    onset: Option<TermId>,
    frequency: Option<Frequency>,
    sex: Option<Sex>,
    modifiers: Vec<TermId>,
    aspect: Option<Aspect>,
    curators: Vec<String>,
    state: PhantomData<(D, P, A)>,
}

impl<P, A> HpoAnnotationBuilder<Unset, P, A> {
    /// Set the disease identifier (e.g. `OMIM:154700`) and its name (e.g. `Marfan syndrome`).
    pub fn disease(
        self,
        disease_id: impl Into<TermId>,
        name: impl ToString,
    ) -> HpoAnnotationBuilder<Set, P, A> {
        HpoAnnotationBuilder {
            disease: Some((disease_id.into(), name.to_string())),
            is_negated: self.is_negated,
            phenotype_term_id: self.phenotype_term_id,
            annotation_references: self.annotation_references,
            onset: self.onset,
            frequency: self.frequency,
            sex: self.sex,
            modifiers: self.modifiers,
            aspect: self.aspect,
            curators: self.curators,
            state: PhantomData,
        }
    }
}

impl<D, A> HpoAnnotationBuilder<D, Unset, A> {
    /// Set the phenotype annotation term id.
    pub fn phenotype(self, phenotype: impl Into<TermId>) -> HpoAnnotationBuilder<D, Set, A> {
        HpoAnnotationBuilder {
            disease: self.disease,
            is_negated: self.is_negated,
            phenotype_term_id: Some(phenotype.into()),
            annotation_references: self.annotation_references,
            onset: self.onset,
            frequency: self.frequency,
            sex: self.sex,
            modifiers: self.modifiers,
            aspect: self.aspect,
            curators: self.curators,
            state: PhantomData,
        }
    }
}

impl<D, P> HpoAnnotationBuilder<D, P, Unset> {
    /// Set the aspect of the HPO annotation.
    pub fn aspect(self, aspect: Aspect) -> HpoAnnotationBuilder<D, P, Set> {
        HpoAnnotationBuilder {
            disease: self.disease,
            is_negated: self.is_negated,
            phenotype_term_id: self.phenotype_term_id,
            annotation_references: self.annotation_references,
            onset: self.onset,
            frequency: self.frequency,
            sex: self.sex,
            modifiers: self.modifiers,
            aspect: Some(aspect),
            curators: self.curators,
            state: PhantomData,
        }
    }

    /// Set aspect to [`Aspect::Phenotype`].
    pub fn aspect_phenotype(self) -> HpoAnnotationBuilder<D, P, Set> {
        self.aspect(Aspect::Phenotype)
    }

    /// Set aspect to [`Aspect::ClinicalCourse`].
    pub fn aspect_clinical_course(self) -> HpoAnnotationBuilder<D, P, Set> {
        self.aspect(Aspect::ClinicalCourse)
    }

    /// Set aspect to [`Aspect::Inheritance`].
    pub fn aspect_inheritance(self) -> HpoAnnotationBuilder<D, P, Set> {
        self.aspect(Aspect::Inheritance)
    }
}

impl<D, P, A> HpoAnnotationBuilder<D, P, A> {
    /// Indicate that the phenotype is *NOT* a characteristic of the annotated disease.
    pub fn negated(mut self) -> Self {
        self.is_negated = true;
        self
    }

    /// Indicate that the phenotype is a characteristic of the annotated disease.
    ///
    /// The phenotype is observed by default and this method does not have to be set.
    pub fn not_negated(mut self) -> Self {
        self.is_negated = false;
        self
    }

    /// Add an annotation reference to the annotation.
    pub fn push_annotation_reference(mut self, annotation_reference: AnnotationReference) -> Self {
        self.annotation_references.push(annotation_reference);
        self
    }

    /// Add multiple annotation references to the annotation.
    pub fn extend_annotation_references(
        mut self,
        annotation_references: impl IntoIterator<Item = AnnotationReference>,
    ) -> Self {
        self.annotation_references.extend(annotation_references);
        self
    }

    /// Clear all previously added annotation references.
    pub fn clear_annotation_references(mut self) -> Self {
        self.annotation_references.clear();
        self
    }

    /// Set the onset of the annotation. The `onset` should be a member of HPO's
    /// [Onset](https://hpo.jax.org/browse/term/HP:0003674) module.
    pub fn onset(mut self, onset: impl Into<TermId>) -> Self {
        self.onset = Some(onset.into());
        self
    }

    /// Set the frequency of the annotation.
    pub fn frequency(mut self, frequency: impl Into<Frequency>) -> Self {
        self.frequency = Some(frequency.into());
        self
    }

    /// Set sex of the annotation.
    pub fn sex(mut self, sex: Sex) -> Self {
        self.sex = Some(sex);
        self
    }

    /// Set sex to [`Sex::Male`].
    pub fn sex_male(self) -> Self {
        self.sex(Sex::Male)
    }

    /// Set sex to [`Sex::Female`].
    pub fn sex_female(self) -> Self {
        self.sex(Sex::Female)
    }

    /// Set sex to [`Sex::Unknown`].
    pub fn sex_unknown(self) -> Self {
        self.sex(Sex::Unknown)
    }

    /// Add a clinical modifier.
    pub fn push_modifier(mut self, modifier: impl Into<TermId>) -> Self {
        self.modifiers.push(modifier.into());
        self
    }

    /// Add several clinical modifiers.
    pub fn extend_modifiers(
        mut self,
        modifiers: impl IntoIterator<Item = impl Into<TermId>>,
    ) -> Self {
        self.modifiers
            .extend(modifiers.into_iter().map(|v| v.into()));
        self
    }

    /// Remove all previously added clinical modifiers.
    pub fn clear_modifiers(mut self) -> Self {
        self.modifiers.clear();
        self
    }

    /// Add a curator record.
    pub fn push_curator(mut self, curator: impl ToString) -> Self {
        self.curators.push(curator.to_string());
        self
    }

    /// Add several curators.
    pub fn extend_curators(mut self, curators: impl IntoIterator<Item = impl ToString>) -> Self {
        self.curators
            .extend(curators.into_iter().map(|v| v.to_string()));
        self
    }

    /// Clear previously added curators.
    pub fn clear_curators(mut self) -> Self {
        self.curators.clear();
        self
    }
}

impl HpoAnnotationBuilder<Set, Set, Set> {
    /// Build the annotation.
    ///
    /// The build is infallible due to static checking in the builder.
    pub fn build(mut self) -> HpoAnnotation {
        let (disease_id, disease_name) = self
            .disease
            .expect("build is callable only after both disease ID and name are set to the builder");

        self.annotation_references.sort_by(|l, r| {
            l.term_id
                .cmp(&r.term_id)
                .then(l.evidence_code.cmp(&r.evidence_code))
        });
        self.modifiers.sort();
        self.curators.sort();

        HpoAnnotation {
            disease_id,
            disease_name,
            is_negated: self.is_negated,
            phenotype_term_id: self
                .phenotype_term_id
                .expect("Build is callable only after phenotype ID is set to the builder"),
            annotation_references: self.annotation_references,
            onset: self.onset,
            frequency: self.frequency,
            sex: self.sex,
            modifiers: self.modifiers,
            aspect: self
                .aspect
                .expect("Build is callable only after aspect is set to the builder"),
            curators: self.curators,
        }
    }
}

/// The HPO annotation corpus with entries parsed to the highest level possible
/// while making no wild assumptions about the parsed data.
#[derive(Debug, Clone, PartialEq)]
pub struct HpoAnnotations {
    /// The HPO annotation records.
    pub lines: Vec<HpoAnnotation>,
    /// The HPOA version (e.g. `2023-04-05`)
    pub version: String,
    pub hpo_version: String,
}

/// Parse disease-phenotype annotations from HPO annotation file.
pub mod io {
    use super::{
        AnnotationReference, Aspect, EvidenceCode, FrequencyParseError, HpoAnnotation,
        HpoAnnotations, Sex,
    };
    use crate::io::{AnnotationWriter, ReadAnnotation};
    use crate::{
        format::Hpoa,
        io::{AnnotationLoadError, AnnotationLoader, ValidationIssue, WriteAnnotation},
    };
    use ontolius::{Prefix, TermIdParseError};
    use regex::Regex;
    use std::collections::BTreeMap;
    use std::fmt::{Debug, Display};
    use std::io::{BufRead, Write};
    use std::sync::LazyLock;

    const HPOA_HEADER: [&str; 12] = [
        "database_id",
        "disease_name",
        "qualifier",
        "hpo_id",
        "reference",
        "evidence",
        "onset",
        "frequency",
        "sex",
        "modifier",
        "aspect",
        "biocuration",
    ];

    // HPOA is a tab-delimited table ...
    const HPOA_DELIMITER: &str = "\t";
    // ... with 12 columns.
    const HPOA_COLUMN_COUNT: usize = HPOA_HEADER.len();
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

    impl Display for HpoaError<'_> {
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
                    FrequencyParseError::PercentageOutOfBounds => write!(
                        f,
                        "Frequency out of bounds: {}",
                        self.fields[FREQUENCY_COL_IDX]
                    ),
                    FrequencyParseError::UnparsableValue => write!(
                        f,
                        "Unparsable frequency value: {}",
                        self.fields[FREQUENCY_COL_IDX]
                    ),
                    FrequencyParseError::NGreaterThanM => write!(
                        f,
                        "`n` is greater than `m`: {}",
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

    /// The reasons for failure in parsing of HPO annotation line in [`HpoAnnotation::read_str`].
    #[derive(Debug, thiserror::Error, PartialEq)]
    pub enum HpoaErrorReason {
        #[error("Cannot parse a record with {0}!={col_cnt}", col_cnt=HPOA_COLUMN_COUNT)]
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

    impl AnnotationWriter<Hpoa> for &HpoAnnotations {
        type Err = std::io::Error;

        fn store_to_write<W>(self, mut write: W) -> Result<(), Self::Err>
        where
            W: Write,
        {
            // Comments
            // #description: "HPO annotations for rare diseases [8362: OMIM; 47: DECIPHER; 4283 ORPHANET]"
            let disease_counts = count_diseases(&self.lines);
            write!(
                &mut write,
                "#description: \"HPO annotations for rare diseases"
            )?;
            if let Some((first, rest)) = disease_counts.split_first() {
                write!(&mut write, " [{}: {}", first.1, first.0)?;

                for (prefix, cnt) in rest {
                    write!(&mut write, "; {}: {}", cnt, prefix)?;
                }

                writeln!(&mut write, "]\"")?;
            } else {
                writeln!(&mut write, "\"")?;
            }
            // versions, tracker
            writeln!(&mut write, "#version: {}", self.version)?;
            writeln!(
                &mut write,
                "#tracker: https://github.com/obophenotype/human-phenotype-ontology/issues"
            )?;
            writeln!(
                &mut write,
                "#hpo-version: https://purl.obolibrary.org/obo/hp/releases/{}/hp.json",
                self.hpo_version
            )?;

            // header
            if let Some((first, rest)) = HPOA_HEADER.split_first() {
                write!(&mut write, "{first}")?;
                for val in rest {
                    write!(&mut write, "{}{}", HPOA_DELIMITER, val)?;
                }
            }
            writeln!(&mut write)?;

            // lines
            for line in &self.lines {
                line.write_ann(&mut write)?;
            }

            Ok(())
        }
    }

    fn count_diseases(lines: &[HpoAnnotation]) -> Vec<(Prefix<'_>, u32)> {
        let mut counts = BTreeMap::new();

        for ann in lines {
            *counts.entry(ann.disease_id.prefix()).or_default() += 1;
        }

        counts.into_iter().collect()
    }

    // #version: 2025-05-06
    static HPOA_VERSION_PT: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"^#(date|version): (?<version>[\w-]+)\w?$")
            .expect("Default pattern should be valid")
    });

    // #hpo-version: http://purl.obolibrary.org/obo/hp/releases/2025-05-06/hp.json
    static HPO_VERSION_PT: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"^#hpo-version: (http|https)://[\w./]+/(?<version>\d{4}-\d{2}-\d{2})/hp\.json$")
            .expect("Default pattern should be valid")
    });

    /// Load the phenotype-disease annotations from HPO annotation file.
    ///
    /// # Example
    ///
    /// ```
    /// use oboannotation::io::AnnotationLoader;
    /// use oboannotation::hpo::HpoAnnotations;
    ///
    /// let data = HpoAnnotations::load_from_path("data/phenotype.real-shortlist.hpoa")
    ///              .expect("The example data should be well formatted");
    ///
    /// // Loaded HPO annotations version `2023-04-05` ...
    /// assert_eq!(data.version.as_str(), "2023-04-05");
    ///
    /// // ... generated with HPO version `2023-04-05` ...
    /// assert_eq!(data.hpo_version.as_str(), "2023-04-05");
    ///
    /// // ... consisting of 86 lines.
    /// assert_eq!(data.lines.len(), 86);
    /// ```
    impl AnnotationLoader<Hpoa> for HpoAnnotations {
        fn load_from_buf_read<R>(mut read: R) -> Result<HpoAnnotations, AnnotationLoadError>
        where
            R: BufRead,
        {
            let mut lines = vec![];
            let mut errors = vec![];
            let mut version = None;
            let mut hpo_version = None;

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
                            } else if let Some(caps) = HPOA_VERSION_PT.captures(line.trim()) {
                                version = Some(caps["version"].to_string());
                            } else if let Some(caps) = HPO_VERSION_PT.captures(line.trim()) {
                                hpo_version = Some(caps["version"].to_string());
                            }
                        } else {
                            // Data
                            match HpoAnnotation::read_str(&line) {
                                Ok(hal) => lines.push(hal),
                                Err(e) => errors.push(ValidationIssue::new(line_number, &e)),
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
            } else {
                match (version, hpo_version) {
                    (Some(v), Some(hv)) => Ok(HpoAnnotations {
                        lines,
                        version: v,
                        hpo_version: hv,
                    }),
                    _ => Err(AnnotationLoadError::Error("Missing version".into())),
                }
            }
        }
    }

    impl WriteAnnotation<Hpoa> for HpoAnnotation {
        fn write_ann<W>(&self, w: &mut W) -> std::io::Result<()>
        where
            W: Write,
        {
            // database_id, disease_name
            write!(w, "{}\t{}\t", self.disease_id, self.disease_name)?;

            // qualifier
            if self.is_negated {
                write!(w, "NOT")?
            }
            write!(w, "\t")?;

            // hpo_id
            write!(w, "{}\t", self.phenotype_term_id)?;

            // reference, evidence
            if let Some((last, rest)) = self.annotation_references.split_last() {
                for ar in rest.iter() {
                    write!(w, "{};", ar.term_id)?;
                }
                write!(w, "{}\t", last.term_id)?;
                write_annotation_reference_evidence_code(w, last.evidence_code)?;
            } else {
                write!(w, "\t")?
            }
            write!(w, "\t")?;

            // onset
            if let Some(onset) = &self.onset {
                write!(w, "{}", onset)?;
            }
            write!(w, "\t")?;

            // frequency
            if let Some(frequency) = &self.frequency {
                write!(w, "{frequency}")?
            }
            write!(w, "\t")?;

            // sex
            if let Some(sex) = &self.sex {
                format_sex(w, sex)?;
            }
            write!(w, "\t")?;

            // modifier
            write_semicolon_separated_array(w, &self.modifiers)?;
            write!(w, "\t")?;

            // aspect
            match &self.aspect {
                Aspect::Phenotype => write!(w, "P")?,
                Aspect::Inheritance => write!(w, "I")?,
                Aspect::ClinicalCourse => write!(w, "C")?,
                Aspect::Modifier => write!(w, "M")?,
                Aspect::PastMedicalHistory => write!(w, "H")?,
            }
            write!(w, "\t")?;

            // biocuration
            write_semicolon_separated_array(w, &self.curators)?;
            writeln!(w,)?;

            Ok(())
        }
    }

    /// Parse HPO annotation line into `HpoAnnotation`.
    ///
    /// # Errors
    ///
    /// Parsing fails on malformed HPO annotation line.
    /// See [`HpoaErrorReason`] for the possible causes.
    ///
    impl ReadAnnotation<Hpoa> for HpoAnnotation {
        type Err = HpoaErrorReason;

        fn read_str(val: &str) -> Result<Self, Self::Err> {
            let fields: Vec<_> = val.trim().split(HPOA_DELIMITER).collect();
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
                        FrequencyParseError::PercentageOutOfBounds
                        | FrequencyParseError::UnparsableValue
                        | FrequencyParseError::NGreaterThanM => {
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

    #[cfg(test)]
    mod test_hpo_io {
        use super::{HpoAnnotation, HpoAnnotations, HpoaErrorReason};
        use crate::format::Hpoa;
        use crate::hpo::io::hpoa_examples::{
            make_complex_hpo_annotation, make_hpo_annotations, make_simple_hpo_annotation,
        };
        use crate::hpo::{AnnotationReference, Aspect, EvidenceCode, FrequencyData};
        use crate::io::{AnnotationLoader, AnnotationWriter, ReadAnnotation, WriteAnnotation};
        use ontolius::TermIdParseError;
        use std::io::BufRead;

        #[test]
        fn write_hpo_annotations() -> Result<(), Box<dyn std::error::Error>> {
            let annotations = make_hpo_annotations();

            let mut w = Vec::new();
            annotations.store_to_write(&mut w)?;

            let lines: Vec<_> = w.lines().map(|l| l.unwrap()).collect();
            assert_eq!(
                &lines,
                &[
                    "#description: \"HPO annotations for rare diseases [4: OMIM; 4: ORPHA]\"",
                    "#version: 2025-05-06",
                    "#tracker: https://github.com/obophenotype/human-phenotype-ontology/issues",
                    "#hpo-version: https://purl.obolibrary.org/obo/hp/releases/2025-05-06/hp.json",
                    "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration",
                    "OMIM:154700\tMarfan syndrome\t\tHP:0002616\tPMID:33436942\tPCS\t\t45/58\t\t\tP\tHPO:probinson[2012-04-24];HPO:probinson[2021-04-01]",
                    "OMIM:154700\tMarfan syndrome\t\tHP:0001647\tPMID:33436942\tPCS\t\t1/58\t\t\tP\tHPO:probinson[2021-04-01]",
                    "OMIM:154700\tMarfan syndrome\t\tHP:0000678\tPMID:33436942\tPCS\t\t8/53\t\t\tP\tHPO:probinson[2012-04-24];HPO:probinson[2021-04-01]",
                    "OMIM:154700\tMarfan syndrome\t\tHP:0008138\tPMID:28050285\tPCS\t\t31/146\t\t\tP\tHPO:probinson[2021-05-27]",
                    "ORPHA:79414\tWoolly hair nevus\t\tHP:0002212\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                    "ORPHA:79414\tWoolly hair nevus\t\tHP:0002213\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                    "ORPHA:79414\tWoolly hair nevus\t\tHP:0011365\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                    "ORPHA:79414\tWoolly hair nevus\t\tHP:0040149\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                ]
            );
            Ok(())
        }

        #[test]
        fn roundtrip_hpo_annotations() -> Result<(), Box<dyn std::error::Error>> {
            let expected = make_hpo_annotations();

            let mut w = Vec::new();
            expected.store_to_write(&mut w)?;

            let actual = HpoAnnotations::load_from_buf_read(&mut &w[..])?;

            assert_eq!(actual, expected);
            Ok(())
        }

        #[test]
        fn read_hpoa_str_line_ok() {
            let line = "OMIM:154700\tMarfan syndrome\t\tHP:0001377\tPMID:28050285;PMID:33436942\tPCS\t\t29/199\t\t\tP\tHPO:probinson[2021-05-27];HPO:probinson[2021-04-01]";

            let hpo_line: Result<HpoAnnotation, _> = HpoAnnotation::read_str(line);

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
                hpo_line.frequency.unwrap().data(),
                &FrequencyData::Ratio { n: 29, m: 199 }
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
        fn read_str_bad_disease_id() {
            let line = "OMIM-154700\tMarfan syndrome\t\tHP:0001377\tPMID:28050285;PMID:33436942\tPCS\t\t29/199\t\t\tP\tHPO:probinson[2021-05-27];HPO:probinson[2021-04-01]";
            let hpo_line: Result<_, _> = HpoAnnotation::read_str(line);

            assert!(hpo_line.is_err());
            assert_eq!(
                hpo_line.unwrap_err(),
                HpoaErrorReason::InvalidDiseaseId(TermIdParseError::MissingDelimiter)
            );
        }

        #[test]
        fn test_write_complex_hpo_annotation() -> std::io::Result<()> {
            let a = make_complex_hpo_annotation();

            let mut buf = Vec::new();
            <HpoAnnotation as WriteAnnotation<Hpoa>>::write_ann(&a, &mut buf)?;

            let line = str::from_utf8(&buf).unwrap();
            assert_eq!(
                line,
                "OMIM:303110\tXq21 deletion syndrome\tNOT\tHP:0000365\tPMID:3476958\tPCS\tHP:0003577\t4/8\tMALE\tHP:0012828\tP\tHPO:iea[2009-02-17];HPO:probinson[2021-09-27];HPO:probinson[2021-09-27]\n"
            );
            Ok(())
        }

        #[test]
        fn test_write_simple_hpo_annotation() -> std::io::Result<()> {
            let a = make_simple_hpo_annotation();

            let mut buf = Vec::new();
            <HpoAnnotation as WriteAnnotation<Hpoa>>::write_ann(&a, &mut buf)?;

            let line = str::from_utf8(&buf).unwrap();
            assert_eq!(
                line,
                "OMIM:303110\tXq21 deletion syndrome\t\tHP:0001419\t\t\t\t\t\t\tI\t\n"
            );
            Ok(())
        }
    }

    #[cfg(test)]
    mod hpoa_examples {
        use crate::hpo::io::HpoAnnotations;
        use crate::hpo::{
            AnnotationReference, Aspect, EvidenceCode, Frequency, HpoAnnotation, Sex,
        };
        use crate::io::ReadAnnotation;

        pub(super) fn make_complex_hpo_annotation() -> HpoAnnotation {
            HpoAnnotation {
                disease_id: ("OMIM", "303110").into(),
                disease_name: "Xq21 deletion syndrome".into(),
                is_negated: true,
                phenotype_term_id: ("HP", "0000365").into(), // Hearing impairment
                annotation_references: vec![AnnotationReference::new(
                    ("PMID", "3476958").into(),
                    EvidenceCode::PCS,
                )],
                onset: Some(("HP", "0003577").into()),
                frequency: Some(Frequency::from_ratio(4u8, 8u8).unwrap()),
                sex: Some(Sex::Male),
                modifiers: vec![("HP", "0012828").into()],
                aspect: Aspect::Phenotype,
                curators: vec![
                    "HPO:iea[2009-02-17]".into(),
                    "HPO:probinson[2021-09-27]".into(),
                    "HPO:probinson[2021-09-27]".into(),
                ],
            }
        }

        pub(super) fn make_simple_hpo_annotation() -> HpoAnnotation {
            HpoAnnotation {
                disease_id: ("OMIM", "303110").into(),
                disease_name: "Xq21 deletion syndrome".into(),
                is_negated: false,
                phenotype_term_id: ("HP", "0001419").into(), // X-linked recessive inheritance
                annotation_references: vec![],
                onset: None,
                frequency: None,
                sex: None,
                modifiers: vec![],
                aspect: Aspect::Inheritance,
                curators: vec![],
            }
        }

        pub(super) fn make_hpo_annotations() -> HpoAnnotations {
            let mut lines = Vec::new();
            lines.extend(marfan_annotations());
            lines.extend(wooly_annotations());
            HpoAnnotations {
                lines,
                version: "2025-05-06".to_string(),
                hpo_version: "2025-05-06".to_string(),
            }
        }
        fn marfan_annotations() -> impl Iterator<Item = HpoAnnotation> {
            [
                "OMIM:154700\tMarfan syndrome\t\tHP:0002616\tPMID:33436942\tPCS\t\t45/58\t\t\tP\tHPO:probinson[2012-04-24];HPO:probinson[2021-04-01]",
                "OMIM:154700\tMarfan syndrome\t\tHP:0001647\tPMID:33436942\tPCS\t\t1/58\t\t\tP\tHPO:probinson[2021-04-01]",
                "OMIM:154700\tMarfan syndrome\t\tHP:0000678\tPMID:33436942\tPCS\t\t8/53\t\t\tP\tHPO:probinson[2012-04-24];HPO:probinson[2021-04-01]",
                "OMIM:154700\tMarfan syndrome\t\tHP:0008138\tPMID:28050285\tPCS\t\t31/146\t\t\tP\tHPO:probinson[2021-05-27]",
            ].into_iter()
                .map(|line| HpoAnnotation::read_str(line).unwrap())
        }

        fn wooly_annotations() -> impl Iterator<Item = HpoAnnotation> {
            [
                "ORPHA:79414\tWoolly hair nevus\t\tHP:0002212\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                "ORPHA:79414\tWoolly hair nevus\t\tHP:0002213\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                "ORPHA:79414\tWoolly hair nevus\t\tHP:0011365\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
                "ORPHA:79414\tWoolly hair nevus\t\tHP:0040149\tORPHA:79414\tTAS\t\tHP:0040281\t\t\tP\tORPHA:orphadata[2025-05-06]",
            ].into_iter()
                .map(|line| HpoAnnotation::read_str(line).unwrap())
        }
    }

    fn format_sex<W>(w: &mut W, sex: &Sex) -> std::io::Result<()>
    where
        W: Write,
    {
        match sex {
            Sex::Unknown => write!(w, "UNKNOWN"),
            Sex::Male => write!(w, "MALE"),
            Sex::Female => write!(w, "FEMALE"),
        }
    }

    fn write_annotation_reference_evidence_code<W>(
        w: &mut W,
        code: EvidenceCode,
    ) -> std::io::Result<()>
    where
        W: Write,
    {
        match code {
            EvidenceCode::IEA => write!(w, "IEA"),
            EvidenceCode::TAS => write!(w, "TAS"),
            EvidenceCode::PCS => write!(w, "PCS"),
        }
    }

    fn write_semicolon_separated_array<W, T>(w: &mut W, a: &[T]) -> std::io::Result<()>
    where
        W: Write,
        T: Display,
    {
        if let Some((last, rest)) = a.split_last() {
            for x in rest {
                write!(w, "{};", x)?;
            }
            write!(w, "{last}")
        } else {
            Ok(())
        }
    }
}
