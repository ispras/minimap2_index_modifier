Minimap2_index_modifier
=======================
Minimap2_index_modifier is a fork of alignment tool [Minimap2](https://github.com/lh3/minimap2).
Unlike the original tool, this can use the variants defined in the VCF file when generating the index, for more accurate alignment.


Minimap2_index_modifier can be used in the same way as the original minimap2. To create a modified index use additional parameter `--vcf-file-with-variants <vcf-file>`.
```bash
minimap2 -d index.mmi --vcf-file-with-variants input.vcf.gz reference.fasta
```

Use flag `--parse-haplotype` if your VCF contains phased haplotypes.

## Contents
* [Installation](#installation)
  * [Compiling from source](#compiling-from-source)
  * [Docker](#docker)
* [Pre-built indexes](#pre-built-indexes)
* [Tests](#tests)

## Installation
### Compiling from source
To compile from source, use this version of tools:

* GCC/G++ 11.4.0+
* HTSlib v1.17

Command to compile:
```bash
cd minimap2_index_modifier && make
```

### Docker
Clone this repository and build a Docker image as follows.
```bash
docker build -t minimap2_index_modifier:2.24 .
```

## Pre-built indexes
This [link](https://nextcloud.ispras.ru/index.php/s/wcb9PpZyr8Gb5CC) contains pre-built modified indexes for next references:
* GRCh38 [(GCA_000001405.15)](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/)
* GRCh37 [(hs37d5)](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh37/)

## Tests
See [test/tests.md](test/tests.md) for more details.

