const _FPATH_SMALL_HPO: &str = "data/hp.small.json";
const FPATH_SMALL_HPOA: &str = "data/phenotype.real-shortlist.hpoa";

mod hpo_annotation_parser {
    use super::FPATH_SMALL_HPOA;
    use oboannotation::hpo::FrequencyData;
    use oboannotation::io::AnnotationWriter;
    use oboannotation::{hpo::HpoAnnotations, io::AnnotationLoader};

    #[test]
    fn load_from_path() {
        let data = HpoAnnotations::load_from_path(FPATH_SMALL_HPOA)
            .expect("Sample data should be well formatted");

        assert_eq!(data.version(), "2023-04-05");
        assert_eq!(data.hpo_version(), "2023-04-05");
        assert_eq!(data.annotations().len(), 86);

        let first = data
            .annotations()
            .first()
            .expect("We should have more than one line");

        assert_eq!(first.disease_id.to_string().as_str(), "OMIM:154700");
        assert_eq!(first.disease_name.as_str(), "Marfan syndrome");
        assert_eq!(first.annotation_references.len(), 2);

        let frequency_data = first.frequency.as_ref().unwrap().data();
        assert_eq!(frequency_data, &FrequencyData::Ratio { n: 29, m: 199 })
    }

    #[test]
    fn write_hpoa_annotations() {
        // Let's assume we obtain the annotations from somewhere.
        let anns = HpoAnnotations::load_from_path(FPATH_SMALL_HPOA)
            .expect("Sample data should be well formatted");

        // We can write them in the HPOA format into a file ...
        // let mut write = File::create("phenotype.new.hpoa").expect("We are allowed to write");
        // ... or into a buffer (this case).
        let mut write = vec![];
        anns.store_to_write(&mut write)
            .expect("Writing into a Vec should not fail");

        let hpoa_lines: Vec<_> = std::str::from_utf8(&write)
            .expect("HPOA should be a utf8 string")
            .lines()
            .collect();

        // The first 5 HPOA lines should look roughly like this ...
        assert_eq!(
            &hpoa_lines[..5],
            &[
                "#description: \"HPO annotations for rare diseases [2: OMIM]\"",
                "#version: 2023-04-05",
                "#tracker: https://github.com/obophenotype/human-phenotype-ontology/issues",
                "#hpo-version: https://purl.obolibrary.org/obo/hp/releases/2023-04-05/hp.json",
                "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration"
            ]
        );

        // ... while the records look like:
        assert_eq!(
            &hpoa_lines[5..8],
            &[
                "OMIM:154700\tMarfan syndrome\t\tHP:0001377\tPMID:28050285;PMID:33436942\tPCS\t\t29/199\t\t\tP\tHPO:probinson[2021-05-27];HPO:probinson[2021-04-01]",
                "OMIM:154700\tMarfan syndrome\t\tHP:0000486\tPMID:8172269\tPCS\t\t110/573\t\t\tP\tHPO:skoehler[2015-07-26];HPO:probinson[2020-08-03]",
                "OMIM:154700\tMarfan syndrome\t\tHP:0005136\tOMIM:154700\tIEA\t\t\t\t\tP\tHPO:probinson[2012-04-24]"
            ]
        );
    }
}

mod create_annotations {
    use oboannotation::hpo::{
        AnnotationReference, Aspect, EvidenceCode, Frequency, FrequencyData, HpoAnnotation,
    };

    #[test]
    fn create_hpoa_entries() {
        // OMIM:154700 Marfan syndrome HP:0001377 PMID:28050285;PMID:33436942 PCS 29/199 P HPO:probinson[2021-05-27];HPO:probinson[2021-04-01]
        let ann = HpoAnnotation::builder()
            .disease(("OMIM", "154700"), "Marfan syndrome")
            .phenotype(("HP", "0001377")) // Limited elbow extension
            .extend_annotation_references([
                AnnotationReference::new(("PMID", "33436942").into(), EvidenceCode::PCS),
                AnnotationReference::new(("PMID", "28050285").into(), EvidenceCode::PCS),
            ])
            .frequency("29/199".parse::<Frequency>().unwrap())
            .aspect_phenotype()
            .extend_curators(["HPO:probinson[2021-05-27]", "HPO:probinson[2021-04-01]"])
            .build();

        assert_eq!(ann.disease_id, ("OMIM", "154700"));
        assert_eq!(ann.disease_name, "Marfan syndrome");
        assert_eq!(ann.is_negated, false);
        assert_eq!(ann.phenotype_term_id, ("HP", "0001377"));
        assert_eq!(ann.annotation_references.len(), 2);
        assert_eq!(
            ann.annotation_references.first().unwrap(),
            &AnnotationReference::new(("PMID", "28050285").into(), EvidenceCode::PCS,)
        );
        assert!(ann.onset.is_none());
        assert_eq!(
            ann.frequency.unwrap().data(),
            &FrequencyData::Ratio { n: 29, m: 199 }
        );
        assert!(ann.sex.is_none());
        assert!(ann.modifiers.is_empty());
        assert_eq!(ann.aspect, Aspect::Phenotype);
        assert_eq!(
            &ann.curators,
            &["HPO:probinson[2021-04-01]", "HPO:probinson[2021-05-27]",]
        );
    }
}
