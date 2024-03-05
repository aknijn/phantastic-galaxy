#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import sys
import os
import argparse
import json
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../PHANtLibs/")
from phantpdf import SampleReport

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phantq_json', dest='phantq_json', help='phantq_json')
    parser.add_argument('--phanta_json', dest='phanta_json', help='phanta_json')
    parser.add_argument('--phantt_json', dest='phantt_json', help='phantt_json')
    parser.add_argument('--phantc_json', dest='phantc_json', help='phantc_json')
    parser.add_argument('--virulotypes', dest='virulotypes', help='virulotypes')
    parser.add_argument('--amrgenes', dest='amrgenes', help='amrgenes')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--strain', dest='strain', help='strain')
    parser.add_argument('--sampledate', dest='sampledate', help='sampledate')
    parser.add_argument('--phantastic_type', dest='phantastic_type', help='phantastic_type')
    parser.add_argument('--samplereport', dest='samplereport', help='samplereport')
    
    args = parser.parse_args()
    report_data = {}
    try:
        phantastic_type = open(args.phantastic_type, 'w')
        # merge JSON files into one
        if os.path.getsize(args.phantq_json) != 0:
            with open(args.phantq_json, "rb") as phantq_infile:
                report_data.update(json.load(phantq_infile))
        if os.path.getsize(args.phanta_json) != 0:
            with open(args.phanta_json, "rb") as phanta_infile:
                report_data.update(json.load(phanta_infile))
        if os.path.getsize(args.phantt_json) != 0:
            with open(args.phantt_json, "rb") as phantt_infile:
                report_data.update(json.load(phantt_infile))
        if os.path.getsize(args.phantc_json) != 0:
            with open(args.phantc_json, "rb") as phantc_infile:
                report_data.update(json.load(phantc_infile))
    finally:
        phantastic_type.write("[" + json.dumps(report_data) + "]")
        phantastic_type.close()

    if len(args.sampledate) < 10:
        sample_date = args.sampledate
    else:
        sample_date = args.sampledate[8:9] + "/" + args.sampledate[5:6] + "/" + args.sampledate[0:3]
    # create sample report
    if "contaminated" not in report_data["qc_messages"]:
        sampleReport = SampleReport(args.species)
        if args.species == "Escherichia coli":
            metadataRow = [report_data["information_name"],report_data["year"],report_data["serotype_o"],report_data["serotype_h"],report_data['qc_status'],
                       report_data["mlst_ST"],report_data["virulotype_stx1"],report_data["virulotype_stx2"],report_data["shigatoxin_subtype"],
                       report_data["virulotype_eae"],report_data["virulotype_ehxa"],sample_date,report_data["coverage"],args.strain]
        else:
            metadataRow = [report_data["information_name"],report_data["year"],report_data['qc_status'],report_data["mlst_ST"],report_data["mlst_CC"],report_data["mlst_lineage"],
                       report_data["serotype_serogroup"],sample_date,report_data["coverage"],args.strain]
        sampleReport.writePdf(metadataRow, args.amrgenes, args.virulotypes, args.samplereport)
        sampleReport.close()
    else:
        Path(args.samplereport).touch()



if __name__ == "__main__":
    main()
