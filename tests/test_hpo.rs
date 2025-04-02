const _FPATH_SMALL_HPO: &str = "data/hp.small.json";
const FPATH_SMALL_HPOA: &str = "data/phenotype.real-shortlist.hpoa";

mod hpo_annotation_parser {

    use oboannotation::{
        hpo::{Frequency, HpoAnnotationLoader},
        io::AnnotationLoader,
    };

    use super::FPATH_SMALL_HPOA;

    #[test]
    fn load_from_path() {
        let parser = HpoAnnotationLoader::default();

        let data = parser
            .load_from_path(FPATH_SMALL_HPOA)
            .expect("Sample data should be well formatted");

        assert_eq!(data.version.as_str(), "2023-04-05");
        assert_eq!(data.lines.len(), 86);

        let first = data
            .lines
            .first()
            .expect("We should have more than one line");

        assert_eq!(first.disease_id.to_string().as_str(), "OMIM:154700");
        assert_eq!(first.disease_name.as_str(), "Marfan syndrome");
        assert_eq!(first.annotation_references.len(), 2);
        assert_eq!(
            first.frequency.as_ref().expect("Should be present"),
            &Frequency::Ratio {
                numerator: 29,
                denominator: 199
            }
        )
    }
}
