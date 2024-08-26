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
import shutil
import subprocess
import datetime
import fileinput
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../PHANtLibs/")
from phantdb import IridaDb
from phantpdf import SampleReport

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def getIdFile(filename):
    splitFilename = filename.split("/")
    if (splitFilename[5][0]=='A'):
         inIdFile = splitFilename[6]
    else:
         inIdFile = splitFilename[5]
    return inIdFile

def getMetadata(inputfiles, inuser, inspecies):
    iridaDb = IridaDb(inspecies)
    with open(inputfiles, 'r') as f:
        content = f.readlines()
    idfiles = [getIdFile(x.rstrip('\n')) for x in content]
    files_id = ",".join(idfiles)
    records = iridaDb.metadata_for_report(inuser, files_id)
    iridaDb.close()
    return records

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='phants filenames')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--phantr_list', dest='phantr_list', help='phantr reports list')
    parser.add_argument('--phantr_reports', dest='phantr_reports', help='phantr reports zip file')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../phantastic.conf')
    IRIDA_DIR = config['fs']['output_path']

    sampleReport = SampleReport(args.species)
    metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), args.species)

    if not os.path.exists('reports'):
        os.mkdir('reports')

    report_list = open(args.phantr_list, 'w')
    report_list.write("Il file reports.zip puo' essere scaricato mediante il pulsante con i tre punti accanto a 'Scarica Tutti i File'\n")
    report_list.write("Il file reports.zip contiene i rapporti di:\n")
    if args.species != "Escherichia coli" and args.species != "Listeria monocytogenes":
        report_list.write("Non e' ancora previsto la creazione di report per questo patogeno\n")
    numColumns = len(sampleReport.dataSommarioHeader)
    for metadataRow in metadata:
        report_list.write(metadataRow[0] + "\n")
        sampleReport.writePdf(metadataRow, IRIDA_DIR + "/" + metadataRow[numColumns + 4], IRIDA_DIR + "/" + metadataRow[numColumns + 3], 'reports/report_' + metadataRow[0] + '.pdf')
    report_list.close()
    shutil.make_archive('irida_reports', format='zip', root_dir='reports')
    shutil.copyfile('irida_reports.zip', args.phantr_reports)
    sampleReport.close()

if __name__ == "__main__":
    __main__()


