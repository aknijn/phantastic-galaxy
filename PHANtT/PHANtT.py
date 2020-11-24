# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import sys
import os
import subprocess
import fileinput

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', dest='species', help='sample species')
    parser.add_argument('--inputs_samples', dest='inputs_samples', help='samples')
    parser.add_argument('--multivirulotypes', dest='multivirulotypes', help='Multi Virulotyping Mapping reads')
    args = parser.parse_args()

    if args.species == "Shiga toxin-producing Escherichia coli":
        thisSpecies = "Escherichia coli"
    else:
        thisSpecies = args.species
    samples_name = []
    samples_virulotypes = []
    # read samples' data from the file in input into three lists
    with open(args.inputs_samples) as f:
        lines = f.read().splitlines()
        for line in lines:
            sample = line.split("\t")
            samples_name.append(sample[0])
            samples_virulotypes.append(sample[1])

    # collect results, individual results in "out_" + sample_name + "_distinct", total genes in genes_found_distinct
    try:
        matrix = open(args.multivirulotypes, 'w')
        if os.path.isfile('genes_found'):
            os.remove('genes_found')
        for (sample_name, sample_virulotypes) in zip(samples_name, samples_virulotypes):
            # select distinct best matching allele & transpose array
            subprocess.call("sort -n -r -k2,2 -k3,3 -k4,4 " + sample_virulotypes  + " | awk '!seen[substr($1, 1, index($1, \"_\")-1)]++ { print $1}' | sort > out_" + sample_name + "_distinct", shell=True)
            # collect all genes found
            subprocess.call("sort -n -r -k2,2 -k3,3 -k4,4 " + sample_virulotypes  + " | awk '!seen[substr($1, 1, index($1, \"_\")-1)]++ { print substr($1, 1, index($1, \"_\")-1)}' >> genes_found", shell=True)
        # filter out duplicate genes
        subprocess.call("awk 'NF && !seen[$1]++ { print $1 }' genes_found | sort > genes_found_distinct", shell=True)
        # create the matrix samples - genes
        with open("genes_found_distinct", "r") as genes_file:
            genes_all = genes_file.read().splitlines()
        matrix.write("sample\t" + "\t".join(genes_all) + "\n")
        for sample_name in samples_name:
            str_sample = sample_name
            with open("out_" + sample_name + "_distinct", "r") as sample_file:
                genes_sample = sample_file.read().splitlines()
                for gene_all in genes_all:
                    this_gene = "-"
                    for gene_sample in genes_sample:
                        if gene_all + "_" in gene_sample:
                            this_gene = gene_sample
                    str_sample = str_sample + "\t" + this_gene
            matrix.write(str_sample + "\n")
    finally:
        matrix.close()
        
    # adjust stx[1/2][a/b] into stx[1/2][A/B] (subunits)
    file = fileinput.FileInput(args.multivirulotypes, inplace=True)
    for line in file:
        line = line.replace('stx1a', 'stx1A')
        line = line.replace('stx2a', 'stx2A')
        line = line.replace('stx1b', 'stx1B')
        line = line.replace('stx2b', 'stx2B')
        print(line, end='')

if __name__ == "__main__":
    __main__()


