## Preface
This is a forked version of the SlamDUNK framework [1]. It was adapted to match our own experimental procedure. 
You can find the original source code [here](https://github.com/t-neumann/slamdunk/releases/latest),
 and read the docs [here](https://t-neumann.github.io/slamdunk/docs.html). 


## Installation for Usage
The code is written in Python. A version 3.6 or above is recommended. It is assumed that `pip` and/or `anaconda`
is installed. Clone this repository and install it via

```commandline
conda env create -f environment.yml
python3 -m pip install .
```

With some few exceptions, the documentation remains the same as already presented
[here](https://t-neumann.github.io/slamdunk/docs.html). However, the major changes include:

1. It is possible to take all reads mapped to all of the genome into account. Generally, this is done by NOT passing
a `bed` file
2. Allow paired mapping. Pass the `--paired` flag to indicate that the two files are paired reads. For example
```commandline
map --paired -5 0 read1.fastq.gz read2.fastq.gz -r reference.fa -o out/
``` 
where `read1.fastq.gz` and `read2.fastq.gz` are your paired reads, `-5 0` means that you take the full 5' end into
account, `-r reference.fa` is your reference file, and `-o out/` is the output directory.

Add `--inverse` to the count or all command if the reads in your library are from the inverse strand. 
### References

[1] Neumann, T., Herzog, V. A., Muhar, M., Haeseler, von, A., Zuber, J., Ameres, S. L., & Rescheneder, P. (2019). [Quantification of experimentally induced nucleotide conversions in high-throughput sequencing datasets](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2849-7). BMC Bioinformatics, 20(1), 258. http://doi.org/10.1186/s12859-019-2849-7


## Original README 
<img src="http://t-neumann.github.io/slamdunk/images/slamdunk_logo_light.png" width="300" title="Slamdunk">

### Streamlining SLAM-Seq analysis with ultra-high sensitivity.

[![GitHub release](https://img.shields.io/github/release/t-neumann/slamdunk.svg)](https://github.com/t-neumann/slamdunk/releases/latest)
[![Travis CI](https://img.shields.io/travis/t-neumann/slamdunk.svg)](https://travis-ci.org/t-neumann/slamdunk)

[![Docker Pulls](https://img.shields.io/docker/pulls/tobneu/slamdunk.svg)](https://hub.docker.com/r/tobneu/slamdunk)
[![Docker Automated build](https://img.shields.io/docker/automated/tobneu/slamdunk.svg)](https://hub.docker.com/r/tobneu/slamdunk/builds/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/slamdunk/README.html)
[![Anaconda build](https://anaconda.org/bioconda/slamdunk/badges/version.svg
)](https://anaconda.org/bioconda/slamdunk)
[![Anaconda downloads](https://anaconda.org/bioconda/slamdunk/badges/downloads.svg
)](https://anaconda.org/bioconda/slamdunk)

[![PyPI release](https://img.shields.io/pypi/v/slamdunk.svg)](https://pypi.python.org/pypi/slamdunk)
![Github Stars](https://img.shields.io/github/stars/t-neumann/slamdunk.svg?style=social&label=Star)

-----

### Slamdunk documentation

http://t-neumann.github.io/slamdunk

### nf-core slamseq workflow

[![nfcore/slamseq](https://github.com/nf-core/slamseq/raw/master/docs/images/nf-core-slamseq_logo.png)](https://nf-co.re/slamseq)

