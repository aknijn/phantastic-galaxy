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

def getMetadata(inputfiles, inuser):
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
    sql = ("select * from v_export_sarscov2 where (email = '" + inuser + "' or right('" + inuser + "',7)='@iss.it') and files_id in (" + files_id + ") order by id")
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

    csv_header = ['id','NomeCampione','DataCampione','DataCaricamento','CodiceInterno','Descrizione','QC','Regione','Provincia','Comune','OrigineIsolato','FasciaEta','NumeroVaccinazioni','DataUltimaVaccinazione','Ospedale','Laboratorio','Lineage','Clade','ORF1ab','Spike','ORF3a','E-protein','M-protein','ORF6','ORF7a','ORF7b','ORF8','N-protein','ORF10','Variante','Sequenza','N_consensus','Coverage']
    args = parser.parse_args()
    metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"))
    if metadata:
        # Create csv
        with open(args.phante_csv, 'w') as meta_csv:
            meta_csv_writer = csv.writer(meta_csv)
            meta_csv_writer.writerow(csv_header)
            for meta_row in metadata:
                meta_csv_writer.writerow(meta_row[0:-2])

if __name__ == "__main__":
    main()
