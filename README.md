# VariantAnnotationTool

A python utility for annotating Variant Call Format (VCF) files with supporting info from the [Broad Institute ExAC Project API](http://exac.hms.harvard.edu/). 

## Dependencies

1. [Python 3.6+](https://www.python.org/downloads/)
2. [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html)

## Install/Setup

1. Install python and confirm version with `python -V`.
2. Install pandas. Using pip: `pip install pandas`.
3. Open terminal in desired install directory and copy/clone git project: `git clone https://github.com/jlbren/vcf-annotation-tool`.

## Usage 
**variant_annotation_tool.py [-h] [--local-only] [-o OUTPUT] input_file**

*positional arguments:*      
* input_file  
  * Input VCF file to be annotated.

*optional arguments:*                  
  * -h, --help  
    * Show this help message and exit
  * --local-only  
    * Optional flag to skip calls to ExAC API.
  * -o OUTPUT, --output  
    * Optional flag to specify output directory. Defaults to current working directory.
## Output

The following table provides a summary of each colomn in the output file. 

| KEY         | DESCRIPTION                                                                              | TYPE   |
|-------------|------------------------------------------------------------------------------------------|--------|
| CHROM       | Chromosome number.                                                                       | int    |
| POS         | Variant nucleotide position.                                                             | int    |
| VARIANT     | Variant allele sequence.                                                                 | string |
| REF         | Reference allele sequence.                                                               | string |
| TYPE        | Mutation type annotation.                                                                | string |
| DPB         | Total read depth per base pair.                                                          | float  |
| AO          | Alternate (variant) allele observation count.                                            | int    |
| AO/RO       | Alternate allele observations / Reference allele observations.                           | float  |
| ALLELE_FREQ | Allele frequency of variant from ExAC Project.                                           | float  |
| GENE        | Effected gene(s) ENSEMBLE ID(s) from ExAC Project. Multiple entries semicolon delimited. | string |
| CONSEQUENCE | Variant consequence effect(s) from ExAC Project. Multiple entries semicolon delimited.   | string |

### Contact

* Jon Brenner | jbrenluc@gmail.com 
