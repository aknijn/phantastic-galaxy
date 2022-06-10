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
import mysql.connector
from mysql.connector import errorcode
from fpdf import FPDF

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))
IRIDA_DIR = '/gfs/irida-phantastic/data/output/'
dataSommarioHeader = ("ID ceppo","Anno","Antigen_O","Antigen_H","QC_status","Regione","MLST_ST","stx1","stx2","stx_subtype","eae","ehxA")
dataVirHeader = ("#gene*", "percentage gene coverage", "gene mean read coverage", "percentage gene identity")
dataAMRHeader = ("Start","Stop","Strand","Gene symbol","Product","Resistance","% Coverage of reference sequence","% Identity to reference sequence")

class PDF(FPDF):
    def header(self):
        # Rendering logos:
        self.image(TOOL_DIR + "/logo_sinistra.png", 20, 8, 26)
        self.image(TOOL_DIR + "/logo_destra.png", 250, 8,26)
        # Setting font
        self.set_font("helvetica", "B", 12)
        # Moving cursor to the right:
        self.cell(120)
        # Printing title:
        self.cell(30, 10, "Laboratorio Nazionale di Riferimento per", align="C")
        self.cell(18)
        self.set_font("helvetica", "BI", 12)
        self.cell(30, 10, "E. coli", align="C")
        # Performing a line break:
        self.ln(10)
        self.cell(120)
        self.set_font("helvetica", "I", 8)
        self.cell(30, 1, "Dipartimento di Sicurezza Alimentare, Nutrizione e Sanità Pubblica Veterinaria", align="C")
        self.ln(4)
        self.cell(120)
        self.cell(30, 1, "Reparto di Sicurezza microbiologica degli alimenti e malattie a trasmissione alimentare", align="C")
        self.ln(6)
        self.cell(120)
        self.set_font("helvetica", "B", 12)
        self.cell(30, 1, "Istituto Superiore di Sanità", align="C")
        self.line(x1=10, y1=36, x2=290, y2=36)
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()} di {{nb}}", align="R")

def getIdFile(filename):
    splitFilename = filename.split("/")
    return splitFilename[5]

def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data

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

    if inspecies == "Escherichia coli":
        sql_species = "SELECT ID_ceppo,Anno,Antigen_O,Antigen_H,QC_status,Regione,MLST_ST,stx1,stx2,stx_subtype,eae,ehxA,DataCampione,Copertura,vir_tab,amr_tab FROM v_report_ecoli"
        sql = (sql_species + " WHERE (email = '" + inuser + "' or right('" + inuser + "',7)='@iss.it') and files_id in (" + files_id + ") group by files_id")
    else:
        sql = "SELECT ID_ceppo,Anno,Antigen_O,Antigen_H,QC_status,Regione,MLST_ST,stx1,stx2,stx_subtype,eae,ehxA,DataCampione,Copertura,vir_tab,amr_tab FROM v_report_ecoli WHERE files_id IN ('0')"

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

def writePdf(dataSommario, dataAMR, dataVir):
    pdf = PDF(orientation="landscape")
    pdf.set_fill_color(r=116, g=183, b=46) 
    pdf.add_page()
    pdf.set_font("helvetica", "B", 14)
    #pdf.cell(10)
    pdf.cell(280, 10 , "Report di sequenziamento genomico di isolati batterici", fill=True, align="C")
    pdf.ln(20)
    pdf.set_font("helvetica", "B", 12)
    pdf.write(8, "Data dell'analisi: " + str(dataSommario[12]))
    pdf.ln(8)
    pdf.write(8, "Copertura del file: " + str(dataSommario[13]) + "X")
    pdf.ln(20)
    pdf.set_font("helvetica", "BU", 14)
    pdf.cell(120)
    pdf.write(8, "Sommario")
    pdf.ln(10)
    pdf.set_font("helvetica", "B", 10)
    line_height = pdf.font_size * 1.5
    col_width = pdf.epw / 11  # distribute content evenly
    pdf.set_fill_color(r=150)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[0], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[1], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[2], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[3], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[4], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    #pdf.multi_cell(col_width, line_height, dataSommarioHeader[5], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[6], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[7], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.set_font("helvetica", "BI", 10)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[8], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[9], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[10], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommarioHeader[11], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)

    pdf.ln(line_height)
    pdf.set_fill_color(r=220)
    pdf.set_font("helvetica", "", 10)
    pdf.multi_cell(col_width, line_height, dataSommario[0], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[1], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[2], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[3], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[4], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    #pdf.multi_cell(col_width, line_height, dataSommario[5], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[6], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[7], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.set_font("helvetica", "I", 10)
    pdf.multi_cell(col_width, line_height, dataSommario[8], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[9], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[10], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataSommario[11], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.ln(line_height)

    pdf.ln(20)
    pdf.set_font("helvetica", "BU", 14)
    pdf.cell(120)
    pdf.write(8, "Virulotipo")
    pdf.ln(10)
    line_height = pdf.font_size * 1.1
    col_width = pdf.epw / 4  # distribute content evenly
    i=1
    pdf.set_fill_color(r=150)
    pdf.set_font("helvetica", "B", 10)
    for cellVirHeader in dataVirHeader:
        pdf.multi_cell(col_width, line_height, cellVirHeader, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, fill=True)
    pdf.ln(line_height)
    for rowVir in dataVir:
        if (i % 2) == 0:
            pdf.set_fill_color(r=255) 
            pdf.set_font("helvetica", "", 10)
        else:
            pdf.set_fill_color(r=220) 
            pdf.set_font("helvetica", "", 10)
        for cellVir in rowVir:
            pdf.multi_cell(col_width, line_height, cellVir, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, fill=True)
        pdf.ln(line_height)
        i = i + 1

    pdf.ln(2)
    pdf.set_font("helvetica", "", 10)
    pdf.write(8, "*=nome del gene_allele_Acc. Number NCBI")
    
    pdf.ln(20)
    pdf.set_font("helvetica", "BU", 14)
    pdf.cell(120)
    pdf.write(8, "Determinanti di antibiotico resistenza")
    pdf.ln(10)
    pdf.set_font("helvetica", "B", 10)
    line_height = pdf.font_size * 1.5
    col_width = pdf.epw / 10  # distribute content evenly
    pdf.set_fill_color(r=150) 
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[0], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[1], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[2], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[3], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(2*col_width, 2*line_height, dataAMRHeader[4], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[5], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, 2*line_height, dataAMRHeader[6], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.multi_cell(col_width, line_height, dataAMRHeader[7], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    pdf.ln(2*line_height)
    for rowAMR in dataAMR:
        if (i % 2) == 0:
            pdf.set_fill_color(r=255) 
            pdf.set_font("helvetica", "", 10)
        else:
            pdf.set_fill_color(r=220) 
            pdf.set_font("helvetica", "", 10)
        for cellAMR in rowAMR:
            pdf.multi_cell(col_width, line_height, cellAMR, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, fill=True)
        pdf.ln(line_height)
        i = i + 1

    # pdf.set_fill_color(r=220)
    # pdf.set_font("helvetica", "", 10)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[0], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[1], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[2], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[3], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(2*col_width, line_height, dataAMR[4], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[5], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[6], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
    # pdf.multi_cell(col_width, 2*line_height, dataAMR[7], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)

    pdf.ln(30)
    pdf.set_font("helvetica", "B", 11)
    pdf.write(8, "Roma " + datetime.date.today().strftime("%d/%m/%Y"))
    pdf.ln(12)
    pdf.cell(30)
    pdf.image(TOOL_DIR + "/firme.png", w=200)
    pdf.output('reports/report_' + dataSommario[0] + '.pdf')

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='phants filenames')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--phantr_reports', dest='phantr_reports', help='phantr reports zip file')
    args = parser.parse_args()
    if args.species == "Shiga toxin-producing Escherichia coli":
        args.species = "Escherichia coli"
    if args.species == "Coronavirus":
        args.species = "SARS-CoV-2"
    metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), args.species)

    if not os.path.exists('reports'):
        os.mkdir('reports')

    for metadataRow in metadata:
        dataSommario = metadataRow[0:13]
        subprocess.call("awk -F '\t' '$4>80 { print $0 }' " + IRIDA_DIR + metadataRow[15] + " | tail -n +2 " + IRIDA_DIR + metadataRow[14] + " > vir_tab_file", shell=True)
        subprocess.call("awk -F '\t' '{ print $3 FS $4 FS $5 FS $6 FS $14 FS $5 FS $10 FS $11 }' " + IRIDA_DIR + metadataRow[15] + " | tail -n +2 > amr_tab_file", shell=True)
        amr_tab = openFileAsTable('amr_tab_file')
        vir_tab = openFileAsTable('vir_tab_file')
        writePdf(dataSommario, amr_tab, vir_tab)
    shutil.make_archive('irida_reports', format='zip', root_dir='reports')
    shutil.copyfile('irida_reports.zip', args.phantr_reports)

if __name__ == "__main__":
    __main__()


