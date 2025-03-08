use std::str::FromStr;

///! Gene product to GO term relations
/// enables links a gene product to a Molecular Function it executes.
/// contributes to links a gene product to a Molecular Function executed by a macromolecular complex, in which the Molecular Function cannot be ascribed to an individual subunit of that complex. Only the complex subunits required for the Molecular Function are annotated to the Molecular Function term with ‘contributes to’.
/// Relations between a gene product and a Biological Process:
/// involved in links a gene product and a Biological Process in which the gene product’s Molecular Function plays an integral role.
/// acts upstream of or within links a gene product and a Biological Process when the mechanism relating the gene product’s activity to the Biological Process is not known.
/// Relations between a gene product and a Cellular Component:
/// is active in links a gene product to the cellular location in which it enables its Molecular Function.
/// located in links a gene product and the Cellular Component, specifically a cellular anatomical anatomy or virion component, in which a gene product has been detected.
/// part of links a gene product and a protein-containing complex.


use serde::Serialize;
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize)]
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
    NegatedAnnotation
}

impl std::fmt::Display for GoTermRelation {
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
            GoTermRelation::NegatedAnnotation => "NOT",
        };
        write!(f, "{}", relation_str)
    }
}

impl FromStr for GoTermRelation {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        if s.starts_with("NOT") {
            return Ok(GoTermRelation::NegatedAnnotation);
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
            _ => Err(format!(
                "Did not recognize '{}' as a GOA relation.",
                s
            )),
        }
    }
}


#[cfg(test)]
mod test {
    use std::assert_eq;

    use super::*;
    // TODO ADD SOME TESTS
}