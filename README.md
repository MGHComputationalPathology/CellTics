# CellTics
Center for Integrated Diagnostics at Mass General Hospital NGS tools

## Installation
```
virtualenv --python=/path/to/python3 venv3
source venv3/bin/activate
python setup.py install
```

## VarGroup
Scans a vcf file and combines multiple nearby SNPs and indels into a single genomic event.

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

### Version 2.0
- added multiprocessing
- converted from python 2.7 to python 3

### Contact

* [Allison MacLeay](mailto:amacleay@mgh.harvard.edu)
* [Ryan Schmidt](mailto:RSCHMIDT@BWH.HARVARD.EDU)
