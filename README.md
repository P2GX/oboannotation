# oboannotation





## Test drive

The following command returns a JSON string with a summary of the GOA file stats (a demo, needs work!).

```shell
cargo run --bin demo -- --goa /../../goa_human.gaf.gz 
```
This should show something like this
```shell
processing /Users/robin/data/go/goa_human.gaf.gz
Parsed 781329 annotations[{"key":"version","value":"24-11-04T20:49"},{"key":"Negated annotations","value":"1494"},{"key":"Total annotations","value":"781329"},{"key":"genes","value":"44541"},{"key":"acts_upstream_of","value":"1061"},{"key":"acts_upstream_of_negative_effect","value":"37"},{"key":"acts_upstream_of_or_within_positive_effect","value":"41"},{"key":"acts_upstream_of_or_within_negative_effect","value":"11"},{"key":"acts_upstream_of_positive_effect","value":"114"},{"key":"colocalizes_with","value":"1016"},{"key":"acts_upstream_of_or_within","value":"2587"},{"key":"part_of","value":"25164"},{"key":"involved_in","value":"189834"},{"key":"located_in","value":"177031"},{"key":"enables","value":"363199"},{"key":"is_active_in","value":"20044"},{"key":"contributes_to","value":"1190"}]
```