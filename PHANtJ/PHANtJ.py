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
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phantq_json', dest='phantq_json', help='phantq_json')
    parser.add_argument('--phanta_json', dest='phanta_json', help='phanta_json')
    parser.add_argument('--phantt_json', dest='phantt_json', help='phantt_json')
    parser.add_argument('--phantc_json', dest='phantc_json', help='phantc_json')
    parser.add_argument('--phantastic_json', dest='phantastic_json', help='phantastic_json')
    
    args = parser.parse_args()
    try:
        report_data = {}
        report = open(args.phantastic_json, 'w')
        # merge JSON files into one
        with open(args.phantq_json, "rb") as phantq_infile:
            report_data.update(json.load(phantq_infile))
        with open(args.phanta_json, "rb") as phanta_infile:
            report_data.update(json.load(phanta_infile))
        with open(args.phantt_json, "rb") as phantt_infile:
            report_data.update(json.load(phantt_infile))
        with open(args.phantc_json, "rb") as phantc_infile:
            report_data.update(json.load(phantc_infile))
    finally:
        report.write("[" + json.dumps(report_data) + "]")
        report.close()

if __name__ == "__main__":
    main()
