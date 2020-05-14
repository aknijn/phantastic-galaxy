#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna

try:
	from allelecall import BBACA
	from utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,uniprot_find,Extract_cgAlleles,RemoveGenes
except ImportError:
	from CHEWBBACA.allelecall import BBACA
	from CHEWBBACA.utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,uniprot_find,Extract_cgAlleles,RemoveGenes

#~ from allelecall import CommonFastaFunctions,callAlleles_protein3,BBACA


def check_if_list_or_folder(folder_or_list):
    list_files = []
    # check if given a list of genomes paths or a folder to create schema
    try:
        f = open(folder_or_list, 'r')
        f.close()
        list_files = folder_or_list
    except IOError:
		
        for gene in os.listdir(folder_or_list):
			
            try:
                genepath = os.path.join(folder_or_list, gene)
                
                if os.path.isdir(genepath):
                    continue
                
                for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print(e)
                pass

    return list_files



def allele_call():

    def msg(name=None):                                                            
        return ''' chewBBACA.py AlleleCall [AlleleCall ...][-h] -i [I] -g [G] -o [O] --cpu [CPU] [-v] [-b [B]][--bsr [BSR]] [--ptf [PTF]] [--fc] [--fr] [--json]
            '''

    parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes when provided a schema",usage=msg())
    #parser.add_argument('AlleleCall', nargs='+', help='do allele call')
    parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
    parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        required=True)
    parser.add_argument("--contained", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("--CDS", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)
    parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False, default=0.6)
    parser.add_argument('--st', nargs='?', type=float, help="size threshold, default at 0.2 means alleles with size variation of +-20 percent will be tagged as ASM/ALM", required=False, default=0.2)
    #parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False, default=False)
    parser.add_argument('--ptf', nargs='?', type=str, help="provide the prodigal training file (ptf) path", required=False, default=False)
    parser.add_argument("--fc", help="force continue", required=False, action="store_true", default=False)
    parser.add_argument("--fr", help="force reset", required=False, action="store_true", default=False)
    parser.add_argument("--json", help="report in json file", required=False, action="store_true", default=False)

    args = parser.parse_args()

    genomeFiles = args.i
    genes = args.g
    cpuToUse = args.cpu
    BSRTresh = args.bsr
    sizeTresh = args.st
    verbose = args.verbose
    BlastpPath = args.b
    gOutFile = args.o
    chosenTaxon = False
    chosenTrainingFile = args.ptf
    forceContinue = args.fc
    forceReset = args.fr
    contained = args.contained
    inputCDS = args.CDS
    jsonReport = args.json

    genes2call = check_if_list_or_folder(genes)

    if isinstance(genes2call, list):
        with open("listGenes2Call.txt", "w") as f:
            for genome in genes2call:
                f.write(genome + "\n")
        genes2call = "listGenes2Call.txt"

    
    #try to open as a fasta
    fasta = SeqIO.parse(genomeFiles, "fasta", generic_dna)
    try:
        isFasta=(any(fasta))
    except:
        isFasta=False
    
    #if is a fasta pass as a list of genomes with a single genome, if not check if is a folder or a txt with a list of paths
    if isFasta==True:
        genomes2call=[os.path.abspath(genomeFiles)]
    else:
        genomes2call = check_if_list_or_folder(genomeFiles)
    
    if isinstance(genomes2call, list):
        with open("listGenomes2Call.txt", "w") as f:
            for genome in genomes2call:
                f.write(genome + "\n")
        genomes2call = "listGenomes2Call.txt"

    
    BBACA.main(genomes2call,genes2call,cpuToUse,gOutFile,BSRTresh,BlastpPath,forceContinue,jsonReport,verbose,forceReset,contained,chosenTaxon,chosenTrainingFile,inputCDS,sizeTresh)


    try:
        os.remove("listGenes2Call.txt")
    except:
        pass
    try:
        os.remove("listGenomes2Call.txt")
    except:
        pass


def main():
    version="2.0.13"
    createdBy="Mickael Silva"
    rep="https://github.com/B-UMMI/chewBBACA"
    contact="mickaelsilva@medicina.ulisboa.pt"

    if len(sys.argv)>1 and "version" in sys.argv[1]:
        print (version)
        return
    
    print ("chewBBACA version "+version+" by "+ createdBy+ " at "+ rep+ "\nemail contact: "+ contact)

    try:
        print ("\n")
        allele_call()
    except Exception as e:
        print (e)
        print('\n\tUSAGE : chewBBACA.py [module] -h \n')
        print('Select one of the following functions :\n')
        i=0
        while i<len(functions_list):
            print (functions_list[i] +" : "+desc_list[i])
            i+=1
if __name__ == "__main__":
    main()
