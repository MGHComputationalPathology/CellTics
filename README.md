# CellTics
Center for Integrated Diagnostics at Mass General Hospital NGS tools

## Installation
```
virtualenv --python=/path/to/python3 venv3
source venv3/bin/activate
python setup.py install
```

If you receive an error pertaining to lzma.h you may need to disable lzma and try python setup.py install again. (This occurs on MacOS Mojave)
```
export HTSLIB_CONFIGURE_OPTIONS=--disable-lzma
python setup.py install
```

With Mac OS High Sierra there is a new security feature to disable multithreading.  If you see an error such as:

```bash
objc[49174]: +[__NSPlaceholderDate initialize] may have been in progress in another thread when fork() was called.                                                                                          
objc[49174]: +[__NSPlaceholderDate initialize] may have been in progress in another thread when fork() was called. We cannot safely call it or ignore it in the fork() child process. Crashing instead. Set 
a breakpoint on objc_initializeAfterForkError to debug.
```

set the following environment variable:
```bash
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

```

## VarGroup
Scans a vcf file and combines multiple nearby SNPs and indels into a single genomic event.

The pickle library DOES NOT work with the cli library which allows calling celltics directly.  When multithreading you must call

```bash
python /path/to/celltics/tools/vargroup.py -i <input> -o <output>
```

VarGrouper runs multithreaded by default.  Use --debug or set threads to 1 (-t 1) to avoid multiprocessing.

### Simplest vargroup command
```bash
celltics vargroup --input-file sorted_variants.vcf --output-file grouped_variants.vcf --ref-seq hg19.fasta -t 1
```

### Run vargroup with bam
```bash
celltics vargroup --input-file sorted_variants.vcf --output-file grouped_variants.vcf --bam-file sorted_alignment.bam --ref-seq hg19.fasta -t 1
```
If a reference sequence is not supplied the UCSC hg19 api is queried ([http://genome.ucsc.edu/](http://genome.ucsc.edu/)).  Variants will be grouped if they are within a certain distance and occur on the same reads.  For more advanced options run ```celltics vargroup --help```.

##Troubleshooting


### Algorithm
![VarGrouper](https://github.com/MGHComputationalPathology/CellTics/blob/master/celltics/docs/graphics/vargrouper_flow.png)

### Version 2.0
- added multiprocessing
- converted from python 2.7 to python 3

### Contact

* [Allison MacLeay](mailto:amacleay@mgh.harvard.edu)
* [Ryan Schmidt](mailto:RSCHMIDT@BWH.HARVARD.EDU)
