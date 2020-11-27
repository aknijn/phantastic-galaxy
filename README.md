# phantastic-galaxy
An automated analysis pipeline PHANtAsTiC: Public Health Analysis of Nucleotides through Assembly, Typing and Clustering

The pipeline is implemented as a Galaxy (https://galaxyproject.org/) workflow integrated with a IRIDA (https://www.irida.ca/) instance.
The tool performs performs trimming (Trimmomatic) and assembly (SPAdes or INNUca), typing and clustering elaborations 
in basis of the type of file (single-end reads versus paired-end reads) and bacterial species (E. coli or Listeria).
In our environment the pipeline is executed on the Galaxy ARIES instance through an API call, 
analysing directly the sequences uploaded in IRIDA.
For more information see: http://biorxiv.org/lookup/doi/10.1101/2020.05.14.095901
