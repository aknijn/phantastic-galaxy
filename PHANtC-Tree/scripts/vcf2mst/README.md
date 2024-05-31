# vcf2mst

Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree.

 - [Link to the Paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0).


When using vcf2mst please use the following citation:

SARS-CoV-2 surveillance in Italy through phylogenomic inferences based on Hamming distances derived from functional annotations of SNPs, MNPs and InDels
Adriano Di Pasquale, Nicolas Radomski, Iolanda Mangone, Paolo Calistri, Alessio Lorusso, Cesare Camma
BMC Genomics 22, 782 (2021). https://doi.org/10.1186/s12864-021-08112-0

## Synopsis

```sh
# Basic usage
vcf2mst.pl input_tsv_file out_file type_of_input [options]

# type_of_input=[vcf,gisaid,algn2pheno,nextclade,tsv,code]

# Examples

# type of inputs
vcf2mst.pl list_of_vcfiles           mst.nwk vcf
vcf2mst.pl gisaid_metadata.tsv       mst.nwk gisaid
vcf2mst.pl algn2pheno_metadata.tsv   mst.nwk algn2pheno
vcf2mst.pl nextclade_metadata.tsv    mst.nwk nextclade
vcf2mst.pl samples_vcfcodes.tsv      mst.nwk code

# profile file only. return a matrix compatible with grapetree input
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -out profile 

# usage 5: filter positions
vcf2mst.pl list_of_vcfiles      mst.tsv     vcf    -minmax-include 10-10000
```

See `examples` folder for the different file formats. See section **Usage** for further details  

# Installation

**Prerequisites**

- Linux OS
- [Perl basic installation](https://www.perl.org/)
- [Graptree](https://github.com/achtman-lab/GrapeTree)

**Installation Step**

- clone/download vcf2mst on your system
- let scripts be executable, e.g., `chmod a+x *.pl`
- put perl scrypts in $PATH, e.g., `cp *.pl /usr/local/bin` 

## Docker based installation

- clone/download vcf2mst on your Linux/Unix based system
- `docker build -t vcf2mst . `


# Usage

```sh
# Basic usage
vcf2mst.pl input_tsv_file out_file type_of_input [options]
```

**type_of_input**=[vcf,gisaid,algn2pheno,nextclade,tsv,code]

## Options

* *-out string=(profile|newick)*: If *string=profile*, the output is just the profile file, without calculating distances and MST. Default is *newick*.
* *-minmax-include value*: Include mutations with position in the "`value`" interval. Format `value` is `min1-max1,min2-max2`. Example `-minmax-include 0-100,200-1500,5000-5500`
* *-minmax-exclude value*: The same as `minmax-include`. mutations in the intervals will be excluded. 
* *-file-minmax-include filename*: Read mutations from filename. Each rows should be a list of intervals with same format of `-minmax-include` value. Comments (with `#`) are allowed in the file.  Example `-file-minmax-include intervals/interval-example.txt`
* *-file-minmax-exclude filename*: The same as `file-minmax-include`. mutations in the intervals will be excluded. 
* *-tsv-XXX*: different options for manipulating a tsv file containing at least 2 columns sample_name and list_of_mutation_codes. This options are considered only in case of *type_of_input=tsv*
  * *-tsv-separator char: the character used as separator on tsv/csv file.(default='\t')
  * *-tsv-sample-pos pos: the position in the tsv file of column containing the sample_name (first position is 0).(default=0)
  * *-tsv-mutationslist-find string=(pos|regexp)*: the way to find the list_of_mutation_codes string in tsv file. if string=pos, -tsv-mutationslist-pos must be set.  if string=regexp, -tsv-mutationslist-regexp must be set. (default=regexp)
  * *-tsv-mutationslist-pos pos*: the position in the tsv file of column containing the list_of_mutation_codes (first position is 0) 
  * *-tsv-mutationslist-regexp string*: the regular expression used to extract the list_of_mutation_codes string. default="`\((.*)\)`"
  * *-tsv-mutation-sep char: the character used as separator between mutations on list_of_mutation_codes string. default=','
  * *-tsv-mutation-pos-regexp string*:  the regular expression used to extract the position of the mutation. default="`^(.*?[_:]?)\w*(\d+)`"
  * *-tsv-mutation-pos-replace string*: the regular expression used to extract the position of the mutation. default="`$1$2`"
* *-debug value*: activate debug with debug_level=*value*. debug is off otherwise (ex: -debug 2)
* *-grapetree-bin string*: string is the grapetree command line. default is `grapetree -p `. An example of docker command is: `docker run  --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p`

**Example of minmax,debug,tsv options**

`perl vcf2mst.pl examples/nextclade_example.tsv profile.tsv tsv -out profile -tsv-sample-pos 0 -tsv-mutationslist-find pos -tsv-mutationslist-pos 15 -tsv-mutation-pos-regexp '\w(\d+)\w' -debug 1 -minmax 9000-10000 -minmax-exclude 9534-9534`

**Example of file-minmax,debug,tsv options**

`perl vcf2mst.pl examples/nextclade_example.tsv profile.tsv tsv -out profile -tsv-sample-pos 0 -tsv-mutationslist-find pos -tsv-mutationslist-pos 15 -tsv-mutation-pos-regexp '\w(\d+)\w'  -file-minmax-include examples/intervals/interval-example-1.txt -file-minmax-exclude examples/intervals/interval-example-2.txt`

**Examples of file intervals format**: see `examples/intervals/interval-example.txt`

## GrapeTree command

In case the vcf2mst is installed locally and the docker version it is not used, there are different way to launch the required `grapetree` tool.

The environment variable GRAPETREE_EXEC can be used for different needs. 
In case there a local installation of grapetree is used, GRAPETREE_EXEC can be set like  `export GRAPETREE_EXEC=grapetree`  (this is the default)
In case there a docker installation of grapetree is used, GRAPETREE_EXEC can be set like  `export GRAPETREE_EXEC=docker run  --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p`  

Grapetree can be also set with `-grapetree-bin` options.

Change the examples according to your specific installation and docker image



# Examples

## From snippy format vcf files 

```sh
# plain
vcf2mst.pl list_of_vcfiles mst.nwk  vcf
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/list_of_vcfiles  /tmp/mst.nwk vcf
``` 

See `examples/list_of_vcfiles` for file format

## From snippy format vcf director

```sh
# plain
vcf2mst.pl folder_with_vcfiles mst.nwk  vcf
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/folder_with_vcfiles /tmp/mst.nwk vcf 
```

See `examples/vcfiles` for folder format 

## From gisaid metadata.tsv 

```sh
# plain
vcf2mst.pl gisaid_metadata.tsv mst.nwk  gisaid
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/gisaid_metadata.tsv /tmp/mst.nwk  gisaid
```

See `examples` folder for  `gisaid_metadata.tsv` file format (minimum information needed) and  `gisaid_full_metadata.tsv` (downloadable from gisaid) file format.


## From vcfcodes csv file

```sh
# plain
vcf2mst.pl samples_vcfcodes.tsv mst.nwk code
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/samples_vcfcodes.tsv  /tmp/mst.nwk code
```

See `examples` folder for `samples_vcfcodes.tsv`  file format

The VCFCODE in samples_vcfcodes.tsv might be:

* a variant code from a vcf file in `snippy` format or
* a GISAID AA_SUBSTITUTIONS value from metadata.tsv (provided by GISAID) or 
* any other variant code produced by any `"variant caller"` of your own choice which maintain the condition `Same variant -> Same code` 

