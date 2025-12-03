const _FPATH_SMALL_HPO: &str = "data/hp.small.json";
const FPATH_SMALL_HPOA: &str = "data/phenotype.real-shortlist.hpoa";

mod hpo_annotation_parser {
    use super::FPATH_SMALL_HPOA;
    use oboannotation::hpo::FrequencyData;
    use oboannotation::{hpo::HpoAnnotations, io::AnnotationLoader};

    #[test]
    fn load_from_path() {
        let data = HpoAnnotations::load_from_path(FPATH_SMALL_HPOA)
            .expect("Sample data should be well formatted");

        assert_eq!(data.version.as_str(), "2023-04-05");
        assert_eq!(data.hpo_version.as_str(), "2023-04-05");
        assert_eq!(data.lines.len(), 86);

        let first = data
            .lines
            .first()
            .expect("We should have more than one line");

        assert_eq!(first.disease_id.to_string().as_str(), "OMIM:154700");
        assert_eq!(first.disease_name.as_str(), "Marfan syndrome");
        assert_eq!(first.annotation_references.len(), 2);

        let frequency_data = first.frequency.as_ref().unwrap().data();
        assert_eq!(frequency_data, &FrequencyData::Ratio { n: 29, m: 199 })
    }
}
