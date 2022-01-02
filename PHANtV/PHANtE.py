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
import mysql.connector
from mysql.connector import errorcode
import shutil
import subprocess
import csv

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))
sql_coronavirus="select sample.id AS id,sample.sampleName AS NomeCampione,date_format(sample.collectionDate,'%Y-%m-%d') AS DataCampione,date_format(sample.arrivalDate,'%Y-%m-%d') AS DataArrivoCampione,date_format(sample.createdDate,'%Y-%m-%d') AS DataCaricamento,sample.strain AS CodiceInterno,sample.description AS Descrizione,max(case when metadata_entry.field_id = 8 then metadata_entry.value end) AS QC,sample.geographicLocationName AS Regione,sample.geographicLocationName2 AS Provincia,sample.geographicLocationName3 AS Comune,sample.isolationSource AS OrigineIsolato,sample.patientAge AS FasciaEta,sample.patientVaccinationNumber AS NumeroVaccinazioni,sample.patientVaccinationDate AS DataUltimaVaccinazione,sample.collectedBy AS Ospedale,user_group.name AS Laboratorio,max(case when metadata_entry.field_id = 7 then metadata_entry.value end) AS Lineage,max(case when metadata_entry.field_id = 54 then metadata_entry.value end) AS Clade,max(case when metadata_entry.field_id = 17 then metadata_entry.value end) AS ORF1ab,max(case when metadata_entry.field_id = 10 then metadata_entry.value end) AS Spike,max(case when metadata_entry.field_id = 9 then metadata_entry.value end) AS ORF3a,max(case when metadata_entry.field_id = 3 then metadata_entry.value end) AS 'E-protein',max(case when metadata_entry.field_id = 16 then metadata_entry.value end) AS 'M-protein',max(case when metadata_entry.field_id = 15 then metadata_entry.value end) AS ORF6,max(case when metadata_entry.field_id = 14 then metadata_entry.value end) AS ORF7a,max(case when metadata_entry.field_id = 13 then metadata_entry.value end) AS ORF7b,max(case when metadata_entry.field_id = 11 then metadata_entry.value end) AS ORF8,max(case when metadata_entry.field_id = 6 then metadata_entry.value end) AS 'N-protein',max(case when metadata_entry.field_id = 18 then metadata_entry.value end) AS ORF10,max(case when metadata_entry.field_id = 35 then metadata_entry.value end) AS Variante,max(case when metadata_entry.field_id = 34 then metadata_entry.value end) AS Sequenza,max(case when metadata_entry.field_id = 52 then if(locate('(',metadata_entry.value) > 0,substr(metadata_entry.value,locate('(',metadata_entry.value) + 1,octet_length(metadata_entry.value) - locate('(',metadata_entry.value) - 2),'') end) AS N_consensus,max(case when qc_entry.total_bases = 0 then 1 else qc_entry.total_bases DIV 29811 end) AS Coverage,lab_users.username AS email,sequence_file_pair_files.files_id AS files_id from metadata_entry join sample on(metadata_entry.sample_id = sample.id) join project_sample on(sample.id = project_sample.sample_id) join user on(user.username = sample.sequencedBy) join user_group_member ugm on(user.id = ugm.user_id) join user_group on(ugm.group_id = user_group.id) join user_group_member lgm on(user_group.id = lgm.group_id) join user lab_users on(lgm.user_id = lab_users.id) join sample_sequencingobject on(sample.id = sample_sequencingobject.sample_id) join sequence_file_pair_files on(sample_sequencingobject.sequencingobject_id = sequence_file_pair_files.pair_id) join qc_entry on(sample_sequencingobject.sequencingobject_id = qc_entry.sequencingObject_id) where project_sample.project_id = 1 and user_group.id > 5 and ugm.role = 'GROUP_MEMBER' and lgm.role = 'GROUP_MEMBER' and "

def getMetadata(inputfiles, inuser, inspecies):
    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../phantastic.conf')
    dbhost = config['db']['host']
    dbdatabase = config['db']['database']
    dbuser = config['db']['user']
    dbpassword = config['db']['password']
    config = {
        'user': dbuser,
        'password': dbpassword,
        'host': dbhost,
        'database': dbdatabase
    }
    with open(inputfiles, 'r') as f:
        content = f.readlines()
    idfiles = [getIdFile(x.rstrip('\n')) for x in content]
    files_id = ",".join(idfiles)
    if inspecies == 'Coronavirus':
        sql = (sql_coronavirus + "(lab_users.username = '" + inuser + "' or right('" + inuser + "',7)='@iss.it') and files_id in (" + files_id + ") group by files_id")
    else:
        sql = "select * from v_export_sarscov2 LIMIT 0"
    try:
        cnx = mysql.connector.connect(**config)
        cursor = cnx.cursor(buffered=True)
        cursor.execute(sql)
        records = cursor.fetchall()
        cursor.close()
        return records
    except mysql.connector.Error as err:
        print(err)
    else:
        cnx.close()

# Obtain idFile from file path
def getIdFile(filename):
    splitFilename = filename.split("/")
    return splitFilename[5]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='input files')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--phante_csv', dest='phante_csv', help='phante_csv')

    csv_header = ['id','NomeCampione','DataCampione','DataArrivoCampione','DataCaricamento','CodiceInterno','Descrizione','QC','Regione','Provincia','Comune','OrigineIsolato','FasciaEta','NumeroVaccinazioni','DataUltimaVaccinazione','Ospedale','Laboratorio','Lineage','Clade','ORF1ab','Spike','ORF3a','E-protein','M-protein','ORF6','ORF7a','ORF7b','ORF8','N-protein','ORF10','Variante','Sequenza','N_consensus','Coverage']
    args = parser.parse_args()
    metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), args.species)
    if metadata:
        # Create csv
        with open(args.phante_csv, 'w') as meta_csv:
            meta_csv_writer = csv.writer(meta_csv)
            meta_csv_writer.writerow(csv_header)
            for meta_row in metadata:
                meta_csv_writer.writerow(meta_row[0:-2])

if __name__ == "__main__":
    main()
