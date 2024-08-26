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
import sys
import os
import shutil
import subprocess
import csv
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../PHANtLibs/")
from phantdb import IridaDb

def getMetadata(inputfiles, inuser, inspecies):
    iridaDb = IridaDb(inspecies)
    with open(inputfiles, 'r') as f:
        content = f.readlines()
    idfiles = [getIdFile(x.rstrip('\n')) for x in content]
    files_id = ",".join(idfiles)
    records = iridaDb.metadata_for_export(inuser, files_id)
    header = iridaDb.header_for_export()
    iridaDb.close()
    return header, records

# Obtain idFile from file path
def getIdFile(filename):
    splitFilename = filename.split("/")
    if (splitFilename[5][0]=='A'):
         inIdFile = splitFilename[6]
    else:
         inIdFile = splitFilename[5]
    return inIdFile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='input files')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--phante_csv', dest='phante_csv', help='phante_csv')

    args = parser.parse_args()
    csv_header, metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), args.species)
    if metadata:
        # Create csv
        with open(args.phante_csv, 'w') as meta_csv:
            meta_csv_writer = csv.writer(meta_csv)
            meta_csv_writer.writerow(csv_header)
            for meta_row in metadata:
                meta_csv_writer.writerow(meta_row[0:-2])

if __name__ == "__main__":
    main()
