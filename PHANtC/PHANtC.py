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
import getopt
import zlib
import textwrap
import shutil


def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', dest='species', help='sample species')
    parser.add_argument('--input', dest='input', help='samples')
    parser.add_argument('--hashprofiles', dest='hashprofiles', help='CRC32 hashes MLST profiles')
    args = parser.parse_args()
    #tool is only kept for compatibility profiles are already CRC32-hashed
    shutil.copyfile(args.input, args.hashprofiles)

if __name__ == "__main__":
    __main__()


