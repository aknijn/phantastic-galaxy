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
import shutil
import subprocess
import datetime
import fileinput
from fpdf import FPDF

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

class SampleReport:
    def __init__(self, species):
        self._species = species
        if species == "Escherichia coli":
            self._dataSommarioHeader = ["ID ceppo ISS","Anno","Antigen O","Antigen H","QC_status","MLST ST","stx1","stx2","stx subtype","eae","ehxA"]
            self._dataSommarioHeaderFormat = ["","","","","","","I","I","I","I","I"]
        elif species == "Listeria monocytogenes":
            self._dataSommarioHeader = ["ID ceppo ISS","Anno","QC status","MLST ST","MLST CC","MLST Lineage","Serogroup"]
            self._dataSommarioHeaderFormat = ["","","","","","",""]
        else:
            self._dataSommarioHeader = ["ID ceppo ISS"]
            self._dataSommarioHeaderFormat = [""]
        self._dataVirHeader = ["#gene*", "percentage gene coverage", "gene mean read coverage", "percentage gene identity"]
        self._dataAMRHeader = ["Start","Stop","Strand","Gene symbol","Product","Resistance","% Coverage of reference sequence","% Identity to reference sequence"]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def species(self):
        return self._species

    @property
    def dataSommarioHeader(self):
        return self._dataSommarioHeader

    def writePdf(self, inSommario, inAMR, inVir, inReportFile):
        subprocess.run("awk -F '\t' '$2>80 { print $0 }' " + inVir + " | tail -n +2 > vir_tab_file", shell=True)
        subprocess.run("sort -k2rn -k4rn -k3rn vir_tab_file | awk '{ if (!seen[substr($1,0,index($1, \"_\"))]++) print $0}' > vir_first_tab_file", shell=True)
        subprocess.run("awk -F '\t' '{ print $3 FS $4 FS $5 FS $6 FS $14 FS $5 FS $10 FS $11 }' " + inAMR + " | tail -n +2 > amr_tab_file", shell=True)
        dataAMR = openFileAsTable('amr_tab_file')
        dataVir = openFileAsTable('vir_first_tab_file')

        numColumns = len(self._dataSommarioHeader)
        if self._species == "Escherichia coli":
            pdf = PDF_STEC(orientation="landscape")
        elif self._species == "Listeria monocytogenes":
            pdf = PDF_LIST(orientation="landscape")
        else:
            pdf = PDF(orientation="landscape")
    
        pdf.set_fill_color(r=116, g=183, b=46) 
        pdf.add_page()
        pdf.set_font("helvetica", "B", 14)
        #pdf.cell(10)
        pdf.cell(280, 10 , "Report di sequenziamento genomico di isolati batterici", fill=True, align="C")
        pdf.ln(20)
        pdf.set_font("helvetica", "B", 12)
        pdf.write(8, "Codice interno del campione: " + str(dataSommario[numColumns + 4]))
        pdf.ln(8)
        pdf.write(8, "Copertura del file: " + str(dataSommario[numColumns + 1]) + "X")
        pdf.ln(8)
        pdf.write(8, "Data dell'analisi: " + str(dataSommario[numColumns]))
        pdf.ln(10)
        pdf.set_font("helvetica", "BU", 14)
        pdf.cell(120)
        pdf.write(8, "Sommario")
        pdf.ln(10)
        pdf.set_font("helvetica", "B", 10)
        line_height = pdf.font_size * 1.5
        col_width = pdf.epw / numColumns  # distribute content evenly
        pdf.set_fill_color(r=150)
        for i in range(len(self._dataSommarioHeader)):
            pdf.set_font("helvetica", "B" + self._dataSommarioHeaderFormat[i], 10)
            pdf.multi_cell(col_width, line_height, self._dataSommarioHeader[i], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
        pdf.ln(line_height)
        pdf.set_fill_color(r=255)
        for i in range(len(self._dataSommarioHeader)):
            pdf.set_font("helvetica", self._dataSommarioHeaderFormat[i], 10)
            pdf.multi_cell(col_width, line_height, dataSommario[i], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
        pdf.ln(line_height)
        pdf.set_font("helvetica", "", 10)
        if len(self._dataSommarioHeader) > 8 and self._dataSommarioHeader[8] == "stx_subtype" and "(*)" in dataSommario[8]:
            pdf.write(8, "(*)=subtype con identità >95% e <100%")
    
        pdf.ln(10)
        pdf.set_font("helvetica", "BU", 14)
        pdf.cell(120)
        pdf.write(8, "Virulotipo")
        pdf.ln(10)
        line_height = pdf.font_size * 1.1
        pdf.cell(25)
        col_width = (pdf.epw - 50)/ 4  # distribute content evenly
        i = 0
        pdf.set_fill_color(r=150)
        pdf.set_font("helvetica", "B", 10)
        for cellVirHeader in self._dataVirHeader:
            pdf.multi_cell(col_width, line_height, cellVirHeader, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, fill=True)
        pdf.ln(line_height)
        for rowVir in dataVir:
            if (i % 2) == 0:
                pdf.set_fill_color(r=255) 
                pdf.set_font("helvetica", "", 10)
            else:
                pdf.set_fill_color(r=220) 
                pdf.set_font("helvetica", "", 10)
            pdf.cell(25)
            j = 0
            for cellVir in rowVir:
                if j == 0:
                    pdf.multi_cell(col_width, line_height, cellVir, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, fill=True)
                else:
                    pdf.multi_cell(col_width, line_height, cellVir, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align="CENTER", fill=True)
                j = 1
            pdf.ln(line_height)
            i += 1
    
        pdf.ln(2)
        pdf.set_font("helvetica", "", 10)
        pdf.cell(25)
        pdf.write(8, "*=nome del gene_allele_Acc. Number NCBI")
    
        pdf.ln(10)
        with pdf.offset_rendering() as dummy:
            # check if title and table header cause page break
            dummy.ln(68)
        if dummy.page_break_triggered:
            # We trigger a page break manually beforehand:
            pdf.add_page()
        pdf.set_font("helvetica", "BU", 14)
        pdf.cell(90)
        pdf.write(8, "Determinanti di antibiotico resistenza")
        pdf.ln(10)
        pdf.set_font("helvetica", "B", 10)
        line_height = pdf.font_size * 1.1
        line_height_list = [3*line_height, 3*line_height, 3*line_height, 3*line_height, 3*line_height, 3*line_height, line_height, line_height]
        col_width = pdf.epw / 8  # distribute columns
        col_width_list = [0.6*col_width, 0.6*col_width, 0.6*col_width, col_width, 2.6*col_width, 0.6*col_width, col_width, col_width]
        col_align_list = ["CENTER", "CENTER", "CENTER", "CENTER", "LEFT", "CENTER", "CENTER", "CENTER"]
        pdf.set_fill_color(r=150) 
        for i in range(8):
            pdf.multi_cell(col_width_list[i], line_height_list[i], self._dataAMRHeader[i], border=1, new_x="RIGHT", new_y="TOP", align="CENTER", fill=True)
        pdf.ln(3*line_height)
        i = 0
        for rowAMR in dataAMR:
            if (i % 2) == 0:
                pdf.set_fill_color(r=255) 
                pdf.set_font("helvetica", "", 10)
            else:
                pdf.set_fill_color(r=220) 
                pdf.set_font("helvetica", "", 10)
            j = 0
            for cellAMR in rowAMR:
                pdf.multi_cell(col_width_list[j], 2*line_height, cellAMR, border=1, new_x="RIGHT", new_y="TOP", max_line_height=pdf.font_size, align=col_align_list[j], fill=True)
                j += 1
            pdf.ln(2*line_height)
            i += 1
    
        pdf.ln(20)
        pdf.set_font("helvetica", "B", 11)
        pdf.write(8, "Roma " + datetime.date.today().strftime("%d/%m/%Y"))
        pdf.ln(12)
        pdf.cell(30)
        if self._species == "Escherichia coli" or self._species == "Listeria monocytogenes":
            pdf.image(TOOL_DIR + "/firme.png", w=200)
        pdf.output(inReportFile)

class PDF(FPDF):
    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()} di {{nb}}", align="R")

class PDF_STEC(FPDF):
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

class PDF_LIST(FPDF):
    def header(self):
        # Rendering logos:
        self.image(TOOL_DIR + "/logo_sinistra.png", 20, 8, 26)
        self.image(TOOL_DIR + "/logo_destra.png", 250, 8,26)
        # Setting font
        self.set_font("helvetica", "B", 12)
        # Moving cursor to the right:
        self.cell(120)
        # Printing title:
        self.cell(30, 10, "Sorveglianza Genomica Listeriosi", align="C")
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

def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data
