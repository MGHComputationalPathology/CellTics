# CellTics
Center for Integrated Diagnostics at Mass General Hospital NGS tools

## Installation
```
virtualenv --python=/path/to/python2.7 venv
source venv/bin/activate
python setup.py install
```

## VarGrouper
A tool to create one vcf row from multiple vcf rows

### Simplest vargroup command
```
celltics vargroup --input-file sorted_variants.vcf --output-file grouped_variants.vcf --ref-seq hg19.fasta
```

### Run vargroup with bam
```
celltics vargroup --input-file sorted_variants.vcf --output-file grouped_variants.vcf --bam-file sorted_alignment.bam --ref-seq hg19.fasta
```
If a reference sequence is not supplied the UCSC hg19 api is queried ([http://genome.ucsc.edu/](http://genome.ucsc.edu/)).  Variants will be grouped if they are within a certain distance and occur on the same reads.  For more advanced options run ```celltics vargroup --help```.

### Algorithm
![VarGrouper](https://github.com/MGHComputationalPathology/CellTics/blob/master/celltics/docs/graphics/vargrouper_flow.png)

### Contact

[Allison MacLeay](mailto:amacleay@mgh.harvard.edu)