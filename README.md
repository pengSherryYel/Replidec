# bayes_classifier
## Aim

use bayes classifier to predict virus lifestyle 

### Important: 

    $PATH should have mmseqs and hmmsearch

## Usage
add the destination dir(where you clone the tool) in your path

```
import sys
sys.path.append("your_path/bayes_classifier")
```

## Example

`from lifestyleBC import bayes_classifier_single,bayes_classifier_batch`

for singe predict    

`bayes_classifier_single("./example/simulate_art_sample1.5.faa","simulate_art_sample1.5","./test")`

output --> list

[prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label]

for batch predict

`bayes_classifier_batch("./example/example.list","./batch_test","BC_predict.summary")`

output --> file same column like above

