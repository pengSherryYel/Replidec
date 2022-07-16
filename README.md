# Replidec: Replication Cycle Detector for Phages
## Aim

Use bayes classifier combine with homology search to predict virus replication cycle

### Important: 

mmseqs, hmmsearch, blastp must set to $PATH, these software can equal or higher than version list below

    MMseqs2 Version: 13.45111

    HMMER 3.3.2 (Nov 2020)

    Protein-Protein BLAST 2.5.0+

```
conda create -n replidec
conda activate replidec
conda install -c bioconda,anaconda "python>=3.8" mmseqs2 hmmer blast
```
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

