#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import configparser
import sys
import os
import csv
import shutil
import subprocess
import json
from md5hash import scan
from pathlib import Path
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../PHANtLibs/")
from phantdb import IridaDb

# Obtain file_id from file path
def get_file_id(filename):
    splitFilename = filename.split("/")
    if (splitFilename[5][0]=='A'):
         inIdFile = splitFilename[6] 
    else:
         inIdFile = splitFilename[5] 
    return inIdFile

def create_allelicprofile_file(allele_database, sample_code, sample_name, phantefsa_allelicprofile):
    allelicprofile = allele_database.allele_strain(sample_code)
    allelicprofile_row = allelicprofile[0]
    allelicprofile_row.replace(sample_code, sample_name, 1)
    with open(phantefsa_allelicprofile, "w") as allelicprofile_file:
        allelicprofile_file.write(allelicprofile_row)

def md5(fname):
    return scan(fname).lower()

def get_experimental_data(inspecies, sample_name, fastq1, fastq2):
    if inspecies == 'Escherichia coli':
        isolate_specie_code = "RF-00003072-PAR"
    elif inspecies == 'Listeria monocytogenes':
        isolate_specie_code = "RF-00000251-MCG"
    elif inspecies == 'Salmonella enterica':
        isolate_specie_code = "RF-00003200-PAR"
    else:
        isolate_specie_code = ""
    if fastq1.endswith('dummy.fastq'):
        library_layout_code = 1
        realfastq1 = os.path.basename(fastq2)
        realfastq1_md5 = md5(fastq2)
    elif fastq2.endswith('dummy.fastq'):
        library_layout_code = 1
        realfastq1 = os.path.basename(fastq1)
        realfastq1_md5 = md5(fastq1)
    else:
        library_layout_code = 2
        realfastq1 = os.path.basename(fastq1)
        realfastq1_md5 = md5(fastq1)
        realfastq2 = os.path.basename(fastq2)
        realfastq2_md5 = md5(fastq2)
    # create json string
    expData = '{'
    expData = ''.join((expData, '"localRawReadId": "', sample_name, '",'))
    expData = ''.join((expData, '"instrumentModelCode": null,'))
    expData = ''.join((expData, '"isolateSpecieCode": "', isolate_specie_code, '",'))
    expData = ''.join((expData, '"serotypeId": null,'))
    expData = ''.join((expData, '"libraryLayoutCode": ', library_layout_code, ','))
    expData = ''.join((expData, '"fastQ1FileName": "', realfastq1, '",'))
    expData = ''.join((expData, '"fastQ1Md5": "', realfastq1_md5, '"'))
    if library_layout_code == 2:
        expData = ''.join((expData, '"fastQ2FileName": "', realfastq2, '"'))
        expData = ''.join((expData, ',"fastQ2Md5": "', realfastq2_md5, '",'))
    expData = ''.join((expData, '}'))
    return expData

def get_ecoli_type_values(phantastic_type):
    with open(phantastic_type, "rb") as phantastic_type_json:
        phantastic_type_dict = json.load(phantastic_type_json)
        coverage = phantastic_type_dict['coverage']
        read_mean_length = phantastic_type_dict['read1_mean_length']
        q30_rate = phantastic_type_dict['q30_rate']
        total_bases = phantastic_type_dict['total_bases']
        serotype = phantastic_type_dict['serotype_o'] + ":" + phantastic_type_dict['serotype_h']
    return (coverage, read_mean_length, q30_rate, total_bases, serotype)

def get_assembly_statistics(phantastic_aq):
    with open(phantastic_aq, "rb") as phantastic_aq_file:
        phantastic_aq_text = phantastic_aq_file.read().splitlines()
    n50_contigs = phantastic_aq_text[10].split("\t")[1]
    genome_size = phantastic_aq_text[7].split("\t")[1]
    number_of_contigs = phantastic_aq_text[5].split("\t")[1]
    return (n50_contigs, genome_size, number_of_contigs)

def get_predicted_pathotype(phantastic_vir):
    # Filter the genes of interest taking care of stx1 and stx2 subtypes
    subprocess.run("awk '{ if (/stx1._|stx2._/ && substr($1, length($1)-1,1)==\"_\" && $2>50) print substr($1, 1, index($1, \"_\")-2) substr($1, length($1),1),$2,$3,$4}' phantastic_vir.tab > virulotype_shortlist")
    subprocess.run("awk '{ if (/stx1._|stx2._/ && substr($1, length($1)-1,1)!=\"_\" && $2>50) print substr($1, 1, index($1, \"_\")-2) \"?\",$2,$3,$4}' phantastic_vir.tab >> virulotype_shortlist")
    subprocess.run("awk '{ if (/eae_|aat_|aggr_|aaic_/ && $2>50) print substr($1, 1, index($1, \"_\")-1),$2,$3,$4}' phantastic_vir.tab >> virulotype_shortlist")
    subprocess.run("sort -n -r -k2,2 -k4,4 -k3,3 virulotype_shortlist | awk !seen[$1]++ {print $1,$2,$4}' > virulotyper_rep", shell=True)
    with open('virulotyper_rep') as virurep:
        virurep.replace("aggr","aggR")
        virurep.replace("aaic","aaiC")
        return virurep

def get_mlst_sequencetype(phantastic_seq):
    with open(phantastic_seq) as phantastic_seq_file:
        phantastic_seq_dict = csv.DictReader(phantastic_seq_file)
    return phantastic_seq_dict

def get_typing_data(inspecies, sample_name, fastq1, fastq2, phantastic_type, phantastic_aq, phantastic_seq, phantastic_vir, predicted_serotype):
    if inspecies == "Escherichia coli":
        (coverage, read_mean_length, q30_rate, total_bases, serotype) = get_ecoli_type_values(phantastic_type)
    (n50_contigs, genome_size, number_of_contigs) = get_assembly_statistics(phantastic_aq)
    virulence_genes = get_predicted_pathotype(phantastic_vir)
    mlst_dict = get_mlst_sequencetype(phantastic_seq)
    
    # create json string
    typeData = '{'
    typeData = ''.join((typeData, '"localRawReadId": "', sample_name, '",'))
    typeData = ''.join((typeData, '"AnalyticalPipelineInfoTag": "PHANtAsTiC V2.1",'))
    typeData = ''.join((typeData, '"QualityCheck": {'))
    typeData = ''.join((typeData, '"Fastp V0.23.2": {'))
    typeData = ''.join((typeData, '"ReadMeanLength": ', read_mean_length, ','))
    typeData = ''.join((typeData, '"Q30Rate": ', q30_rate, ','))
    typeData = ''.join((typeData, '"TotalBases": ', total_bases))
    typeData = ''.join((typeData, '},'))
    typeData = ''.join((typeData, '"AssemblyQualityStatistics": {'))
    typeData = ''.join((typeData, '"AssemblyCoverage": ', + coverage, ','))
    typeData = ''.join((typeData, '"AssemblyStatistics": {'))
    typeData = ''.join((typeData, '"N50Contigs": ', n50_contigs, ','))
    typeData = ''.join((typeData, '"GenomeSize": ', genome_size, ','))
    typeData = ''.join((typeData, '"NumberOfContigs": ', number_of_contigs))
    typeData = ''.join((typeData, '}}},'))
    typeData = ''.join((typeData, '"Results": {'))
    typeData = ''.join((typeData, '"ParamCode": {'))
    typeData = ''.join((typeData, '"PredictedPathotype": {'))
    typeData = ''.join((typeData, '"Pathotype": "NA",'))
    typeData = ''.join((typeData, '"VeroToxin": {'))
    typeData = ''.join((typeData, '"VT1Positive": ', str("stx1" in virulence_genes).lower(), ','))
    typeData = ''.join((typeData, '"VT2Positive": ', str("stx2" in virulence_genes).lower()))
    typeData = ''.join((typeData, '},'))
    typeData = ''.join((typeData, '"AdhesionGenes": {'))
    typeData = ''.join((typeData, '"eaePositive": ', str("eae" in virulence_genes).lower(), ','))
    typeData = ''.join((typeData, '"aatPositive": ', str("aat" in virulence_genes).lower(), ','))
    typeData = ''.join((typeData, '"aggRPositive": ', str("aggR" in virulence_genes).lower(), ','))
    typeData = ''.join((typeData, '"aaiCPositive": ', str("aaiC" in virulence_genes).lower()))
    typeData = ''.join((typeData, ' },'))
    typeData = ''.join((typeData, ' "Software": "patho_typing",'))
    typeData = ''.join((typeData, '"GeneList": ['))
    i = 0
    stx1_subtype = ''
    stx2_subtype = ''
    for row in virulence_genes:
        i += 1
        gene = [str(col).rstrip() for col in row.split('\t')]
        typeData = ''.join((typeData, '{'))
        typeData = ''.join((typeData, '"GeneName": "', gene[0], '",'))
        typeData = ''.join((typeData, '"Identity": ', gene[2], ','))
        typeData = ''.join((typeData, '"Coverage": ', gene[1]))
        if i < len(virulence_genes):
            typeData = ''.join((typeData, '},'))
        else:
            typeData = ''.join((typeData, '}'))
        if "stx1" in gene[0]:
            stx1_subtype = ','.join((stx1_subtype, gene[0]))
        if "stx2" in gene[0]:
            stx2_subtype = ','.join((stx2_subtype, gene[0]))
    if stx1_subtype == '':
        stx1_subtype = 'NT'
    if stx2_subtype == '':
        stx2_subtype = 'NT'
    typeData = ''.join((typeData, ']'))
    typeData = ''.join((typeData, '},'))
    typeData = ''.join((typeData, '"PredictedSerotype": {'))
    typeData = ''.join((typeData, '"Serotype": "', serotype, '",'))
    typeData = ''.join((typeData, '"Software": "PHANtAsTiC V2.1"'))
    typeData = ''.join((typeData, '},'))
    typeData = ''.join((typeData, '"PredictedStxType": "' + ':'.join((stx1_subtype, stx2_subtype)) + '",'))
    typeData = ''.join((typeData, '"MLSTSequenceType": {'))
    typeData = ''.join((typeData, '"ST": ', mlst_dict.get("ST"), ','))
    typeData = ''.join((typeData, '"Software": "mlst V2.16.1",'))
    typeData = ''.join((typeData, '"GeneList": ['))
    mlst_dict.pop("ST")
    typeData = ''.join((typeData, json.dumps(mlst_dict)))
    typeData = ''.join((typeData, ']'))
    typeData = ''.join((typeData, '}}}}'))
    return typeData

def get_epidemiological_data(inspecies, sampleName):
    stecdb = StecDb(args.species)
    metadata = stecdb.metadata_for_efsa(file_ids)
    #sample_name, sampId, sampPoint, sampCountry, origCountry, sampArea, sampY, sampM, sampD, sampling_matrix_code, sampling_matrix_free_text, isolId, YEAR(DateOfSampling), MONTH(DateOfSampling), DAY(DateOfSampling)
    meta_row = metadata[0]
    # create json string
    epiData = '{'
    epiData = ''.join((epiData, '"localRawReadId": "', sample_name, '",'))
    epiData = ''.join((epiData, '"samplingData": {"sampling_local_id": ', meta_row[1], '",'))
    epiData = ''.join((epiData, '"programme_type_id": null,'))
    epiData = ''.join((epiData, '"programme_info_id": null,'))
    epiData = ''.join((epiData, '"sampler_id": "', meta_row[2], '",'))
    epiData = ''.join((epiData, '"sampling_point_id": "', meta_row[3], '",'))
    epiData = ''.join((epiData, '"sampling_country_id": "', meta_row[4], '",'))
    epiData = ''.join((epiData, '"sampling_country_origin_id": "', meta_row[5], '",'))
    epiData = ''.join((epiData, '"sampling_area_id": "', meta_row[6], '",'))
    epiData = ''.join((epiData, '"sampling_year": ', meta_row[7], ','))
    epiData = ''.join((epiData, '"sampling_month": ', meta_row[8], ','))
    epiData = ''.join((epiData, '"sampling_day": ', meta_row[9], ','))
    epiData = ''.join((epiData, '"sampling_matrix_code": "', meta_row[10], '",'))
    epiData = ''.join((epiData, '"sampling_matrix_free_text": "', meta_row[11], '",'))
    epiData = ''.join((epiData, '},'))
    epiData = ''.join((epiData, '"isolationData": {'))
    epiData = ''.join((epiData, '"isolation_local_id": "', meta_row[12], '",'))
    epiData = ''.join((epiData, '"isolation_year": ', meta_row[13], ','))
    epiData = ''.join((epiData, '"isolation_month": ', meta_row[13], ','))
    epiData = ''.join((epiData, '"isolation_day": ', meta_row[13], '}'))
    epiData = ''.join((epiData, '}'))
    return epiData

def create_payload_file(experimental_data, typing_data, epidemiological_data, phantefsa_payload):
    try:
        payload_data = {}
        payload_file = open(phantefsa_payload, 'w')
        # merge JSON files into one
        payload_data.update(json.load(experimental_data))
        payload_data.update(json.load(typing_data))
        payload_data.update(json.load(epidemiological_data))
    finally:
        payload_file.write(json.dumps(payload_data))
        payload_file.close()
    return payload_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='input files')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--phantefsa_payloads', dest='phantefsa_payloads', help='phantefsa_payloads')
    parser.add_argument('--phantefsa_allelicprofiles', dest='phantefsa_allelicprofiles', help='phantefsa_allelicprofiles')
    parser.add_argument('--phantefsa_apiresponses', dest='phantefsa_apiresponses', help='phantefsa_apiresponses')

    args = parser.parse_args()
    if args.species == "Escherichia coli":
        EFSA_ROLE = "EURL STEC"
    else:
        EFSA_ROLE = "XXXX"
    iridadb = IridaDb(args.species)
    if iridadb.user_in_role(args.user.replace("__at__", "@"), EFSA_ROLE):
        with open(args.input_files, 'r') as f:
            content = f.readlines()
        idfiles = [get_file_id(x.rstrip('\n')) for x in content]
        file_ids = ",".join(idfiles)
        metadata = iridadb.metadata_for_efsa(file_ids)
        #sample_id, sample_code, sample_name, fastq1, fastq2, phantastic_type, phantastic_aq, phantastic_seq, phantastic_vir
        try:
            # initialise the cumulative output files
            payload_data = {}
            payloads = open(args.phantefsa_payloads, 'w')
            shutil.copyfile(iridadb.allele_header_file, args.phantefsa_allelicprofiles)
            Path(args.phantefsa_apiresponses).touch()
            if metadata:
                # Upload each sample to the EFSA WGS database
                for meta_row in metadata:
                    # define allelicprofile filename and copy the header row into it
                    allelicprofile = meta_row[2] + ".tsv"
                    shutil.copyfile(iridadb.allele_header_file, allelicprofile)
                    create_allelicprofile_file(iridadb, meta_row[1], meta_row[2], allelicprofile)
                    # save result in the cumulative allelicprofiles file
                    subprocess.run("tail -n 1 " + allelicprofile + " >> " + args.phantefsa_allelicprofiles, shell=True)
                    # obtain the json data for the three parts that will be fused into one payload file
                    payload = meta_row[2] + ".json"
                    experimental_data = get_experimental_data(args.species, meta_row[2], meta_row[3], meta_row[4])
                    typing_data = get_typing_data(args.species, meta_row[3], meta_row[4], meta_row[5], meta_row[6], meta_row[7], meta_row[8])
                    epidemiological_data = get_epidemiological_data(args.species, meta_row[2])
                    phantefsa_payload_json = create_payload_file(experimental_data, typing_data, epidemiological_data, payload)
                    # save payload data for the cumulative json file
                    payload_data.update(json.load(phantefsa_payload_json))
                    # call the EFSA APIs to upload the sample data and save the analyticalPipelineRunId as externalId
                    apiresponse = meta_row[2] + ".txt"
                    analyticalPipelineRunId = call_upload_results(payload, apiresponse)
                    call_upload_results_allelicprofile(allelicprofile, analyticalPipelineRunId, apiresponse)
                    iridadb.update_externalId(analyticalPipelineRunId, meta_row[0])
                    # save response data in the cumulative apiresponses file
                    subprocess.run("echo ================================================= >> " + args.phantefsa_apiresponses, shell=True)
                    subprocess.run("echo " + meta_row[2] + " >> " + args.phantefsa_apiresponses, shell=True)
                    subprocess.run("cat " + apiresponse + " >> " + args.phantefsa_apiresponses, shell=True)
        finally:
            payloads.write("[" + json.dumps(payload_data) + "]")
            payloads.close()
    else:
        subprocess.run("echo 'Utente non autorizzato' > " + args.phantefsa_apiresponses, shell=True)
    iridadb.close()

if __name__ == "__main__":
    main()

