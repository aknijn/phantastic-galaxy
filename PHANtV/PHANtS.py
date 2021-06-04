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
import datetime
import fileinput
import mysql.connector
from mysql.connector import errorcode
import pandas as pd
import numpy as np

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def insertFile(filename, report):
    with open(filename) as html_in:
        for line in html_in:
            report.write(line)

def getIdFile(filename):
    splitFilename = filename.split("/")
    return splitFilename[5]

def getMetadata(inputfiles, species):
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

    if species == "Escherichia coli":
        sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_ecoli WHERE files_id IN (" + files_id + ")")
    elif species == "Listeria monocytogenes":
        sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_listeria WHERE files_id IN (" + files_id + ")")
    elif species == "SARS-CoV-2":
        sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_sarscov2 WHERE files_id IN (" + files_id + ")")
    else:
        sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_ecoli WHERE files_id IN (" + files_id + ")")

    try:
        cnx = mysql.connector.connect(**config)
        cursor = cnx.cursor(buffered=True)
        cursor.execute(sql)
        records = pd.DataFrame(cursor.fetchall(),columns=['Anno','MLST','QC','Regione','Sero','SeroT','Stx12','StxSub','Eae'])
        cursor.close()
        return records
    except mysql.connector.Error as err:
        print(err)
    else:
        cnx.close()

def getPassedFailed(dataframe):
    numPassed = 0
    numFailed = 0
    numPassedFailed = dataframe['QC'].value_counts()
    dfkeys = dataframe['QC'].value_counts().keys().tolist()
    if "Passed" in dfkeys:
        numPassed = numPassedFailed['Passed']
    if "Failed" in dfkeys:
        numFailed = numPassedFailed['Failed']
    strPassedFailed = "var QCData = {labels: ['Passed', 'Failed'],datasets: [{backgroundColor: ['rgb(75, 192, 192)','rgb(255, 159, 64)'],data: [%s,%s]}]};\n" % (numPassed, numFailed)
    return strPassedFailed

def getRegionYear(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Anno"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRY = "var RegioneAnnoData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lstyears = []
    for idx, row in dfpivot.iterrows():
        year = idx
        lstyeardata = []
        for dfkey in dfkeys:
            lstyeardata.append(row[dfkey])
        stryeardata = ",".join([str(item) for item in lstyeardata])
        lstyears.append("{label: \"" + year + "\", backgroundColor : getRCH(), data: [" + stryeardata + "]}")
    strRY = strRY + ",".join(lstyears) + "]};\n"
    return strRY

def getHospitalYear(dataframe, region, regionData):
    regiondataframe = dataframe[dataframe.Regione==region]
    dfpivot = pd.pivot_table(regiondataframe,index=["Anno"], columns='Sero', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strHY = "var " + regionData + " = {labels: [" + ",".join(["\"" + item.replace("\"", "") + "\"" for item in dfkeys]) + "],datasets: ["
    lstyears = []
    for idx, row in dfpivot.iterrows():
        year = idx
        lstyeardata = []
        for dfkey in dfkeys:
            lstyeardata.append(row[dfkey])
        stryeardata = ",".join([str(item) for item in lstyeardata])
        lstyears.append("{label: \"" + year + "\", backgroundColor : getRCH(), data: [" + stryeardata + "]}")
    strHY = strHY + ",".join(lstyears) + "]};\n"
    return strHY

def getST(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], values='QC', aggfunc=len, fill_value=0)
    dfSTs = dfpivot.index
    strST = "var STData = {labels: [" + ",".join(["'" + item + "'" for item in dfSTs]) + "],datasets: [{backgroundColor: ["
    strST = strST + ",".join(["getRCH()" for item in dfSTs]) + "], data: ["
    lstSTdata = []
    for idx, row in dfpivot.iterrows():
        lstSTdata.append(str(row[0]))
    strST = strST + ",".join(lstSTdata) + "]}]};\n"
    return strST

def getSero(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Sero"], values='QC', aggfunc=len, fill_value=0)
    dfSGs = dfpivot.index
    strSG = "var SeroData = {labels: [" + ",".join(["'" + item + "'" for item in dfSGs]) + "],datasets: [{backgroundColor: ["
    strSG = strSG + ",".join(["getRCH()" for item in dfSGs]) + "], data: ["
    lstSGdata = []
    for idx, row in dfpivot.iterrows():
        lstSGdata.append(str(row[0]))
    strSG = strSG + ",".join(lstSGdata) + "]}]};\n"
    return strSG

def getStxSub(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["StxSub"], values='QC', aggfunc=len, fill_value=0)
    dfSXs = dfpivot.index
    strSX = "var StxData = {labels: [" + ",".join(["'" + item + "'" for item in dfSXs]) + "],datasets: [{backgroundColor: ["
    strSX = strSX + ",".join(["getRCH()" for item in dfSXs]) + "], data: ["
    lstSXdata = []
    for idx, row in dfpivot.iterrows():
        lstSXdata.append(str(row[0]))
    strSX = strSX + ",".join(lstSXdata) + "]}]};\n"
    return strSX

def getSeroYear(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Anno"], columns='Sero', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strSY = "var SeroAnnoData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lstyears = []
    for idx, row in dfpivot.iterrows():
        year = idx
        lstyeardata = []
        for dfkey in dfkeys:
            lstyeardata.append(row[dfkey])
        stryeardata = ",".join([str(item) for item in lstyeardata])
        lstyears.append("{label: \"" + year + "\", backgroundColor : getRCH(), data: [" + stryeardata + "]}")
    strSY = strSY + ",".join(lstyears) + "]};\n"
    return strSY

def getRegionSero(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Sero"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRS = "var RegioneSeroData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lstseros = []
    for idx, row in dfpivot.iterrows():
        sero = idx
        lstserodata = []
        for dfkey in dfkeys:
            lstserodata.append(row[dfkey])
        strserodata = ",".join([str(item) for item in lstserodata])
        lstseros.append("{label: \"" + sero + "\", backgroundColor : getRCH(), data: [" + strserodata + "]}")
    strRS = strRS + ",".join(lstseros) + "]};\n"
    return strRS

def getRegionST(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRS = "var RegioneSTData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lststs = []
    for idx, row in dfpivot.iterrows():
        st = idx
        lststdata = []
        for dfkey in dfkeys:
            lststdata.append(row[dfkey])
        strstdata = ",".join([str(item) for item in lststdata])
        lststs.append("{label: \"" + st + "\", backgroundColor : getRCH(), data: [" + strstdata + "]}")
    strRS = strRS + ",".join(lststs) + "]};\n"
    return strRS

def getRegionStx12(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Stx12"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRStx = "var RegioneStx12Data = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lststxs = []
    for idx, row in dfpivot.iterrows():
        stx = idx
        lststxdata = []
        for dfkey in dfkeys:
            lststxdata.append(row[dfkey])
        strstxdata = ",".join([str(item) for item in lststxdata])
        lststxs.append("{label: \"" + stx + "\", backgroundColor : getRCH(), data: [" + strstxdata + "]}")
    strRStx = strRStx + ",".join(lststxs) + "]};\n"
    return strRStx

def getStx12Sero(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Sero"], columns='Stx12', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strSS = "var StxSeroData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lstseros = []
    for idx, row in dfpivot.iterrows():
        sero = idx
        lstserodata = []
        for dfkey in dfkeys:
            lstserodata.append(row[dfkey])
        strserodata = ",".join([str(item) for item in lstserodata])
        lstseros.append("{label: \"" + sero + "\", backgroundColor : getRCH(), data: [" + strserodata + "]}")
    strSS = strSS + ",".join(lstseros) + "]};\n"
    return strSS
    
def getEaeSero(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Sero"], columns='Eae', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strES = "var EaeSeroData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lstseros = []
    for idx, row in dfpivot.iterrows():
        sero = idx
        lstserodata = []
        for dfkey in dfkeys:
            lstserodata.append(row[dfkey])
        strserodata = ",".join([str(item) for item in lstserodata])
        lstseros.append("{label: \"" + sero + "\", backgroundColor : getRCH(), data: [" + strserodata + "]}")
    strES = strES + ",".join(lstseros) + "]};\n"
    return strES

def getRegionLineage(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRL = "var RegioneLineageData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lstlineages = []
    for idx, row in dfpivot.iterrows():
        lineage = idx
        lstlineagedata = []
        for dfkey in dfkeys:
            lstlineagedata.append(row[dfkey])
        strlineagedata = ",".join([str(item) for item in lstlineagedata])
        lstlineages.append("{label: \"" + lineage + "\", backgroundColor : getRCH(), data: [" + strlineagedata + "]}")
    strRL = strRL + ",".join(lstlineages) + "]};\n"
    return strRL

def getRegionMutations(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["SeroT"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strRM = "var RegioneMutazioniData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lstlineages = []
    for idx, row in dfpivot.iterrows():
        mutazioni = idx
        lstlineagedata = []
        for dfkey in dfkeys:
            lstlineagedata.append(row[dfkey])
        strlineagedata = ",".join([str(item) for item in lstlineagedata])
        lstlineages.append("{label: \"" + mutazioni + "\", backgroundColor : getRCH(), data: [" + strlineagedata + "]}")
    strRM = strRM + ",".join(lstlineages) + "]};\n"
    return strRM

def getSTYear(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], columns='Anno', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    strSY = "var STAnnoData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lstlineages = []
    for idx, row in dfpivot.iterrows():
        lineages = idx
        lstlineagedata = []
        for dfkey in dfkeys:
            lstlineagedata.append(row[dfkey])
        strlineagedata = ",".join([str(item) for item in lstlineagedata])
        lstlineages.append("{label: \"" + lineages + "\", backgroundColor : getRCH(), data: [" + strlineagedata + "]}")
    strSY = strSY + ",".join(lstlineages) + "]};\n"
    return strSY

def getYear(dataframe):
    lstyrdata = dataframe.Anno.unique()
    lstyrdata.sort()
    if dataframe.Anno.nunique(dropna = True) == 1:
        strYR = "Anno " + str(lstyrdata[0])
    else:
        strYR = "Anni " + ", ".join([str(item) for item in lstyrdata])
    return strYR

def getTable(dataframe):
    df_clusters = dataframe.groupby('Regione')
    #Column Count
    df_count = df_clusters.count()[['QC']]
    df_count.columns = ['Numero di sequenziamenti']
    strCT = df_count.to_html(classes='table table-cross', border=0)
    return strCT

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='phants filenames')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--phants_stat', dest='phants_stat', help='phants stat report html file')
    parser.add_argument('--phants_trend', dest='phants_trend', help='phants trend report csv file')
    args = parser.parse_args()
    if args.species == "Shiga toxin-producing Escherichia coli":
        args.species = "Escherichia coli"
    if args.species == "Coronavirus":
        args.species = "SARS-CoV-2"
    metadata = getMetadata(args.input_files, args.species)
    # write csv
    dfpivot = pd.pivot_table(metadata,index=["MLST"], columns='Anno', values='QC', aggfunc=len, fill_value=0)
    if args.species == "SARS-CoV-2":
        dfpivot = dfpivot.rename_axis("Lineage")
    dfpivot.to_csv(args.phants_trend) 
    # write html
    try:
        report = open(args.phants_stat, 'w')
        # write head html
        if args.species == "Escherichia coli":
            insertFile(TOOL_DIR + "/report_head_stec.html", report)
        else:
            if args.species == "SARS-CoV-2":
                insertFile(TOOL_DIR + "/report_head_sc2.html", report)
                report.write(getTable(metadata))
                report.write("</td></tr></table>\n")
            else:
                insertFile(TOOL_DIR + "/report_head.html", report)
        if args.species == "Escherichia coli":
            report.write(getYear(metadata))
            if metadata.Anno.nunique(dropna = True) == 1:
                insertFile(TOOL_DIR + "/report_head2_stec1.html", report)
            else:
                insertFile(TOOL_DIR + "/report_head2_stec2.html", report)
        report.write("<div>ISS: Riepilogo per <i>%s</i>, elaborato %s</div>\n<script>\n" % (args.species, datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")))
        # Add the data
        report.write(getPassedFailed(metadata))
        report.write(getRegionYear(metadata))
        report.write(getST(metadata))
        if args.species == "SARS-CoV-2":
            report.write(getRegionLineage(metadata))
            report.write(getRegionMutations(metadata))
            report.write(getRegionStx12(metadata))
            report.write(getHospitalYear(metadata,'Abruzzo','AbData'))
            report.write(getHospitalYear(metadata,'Basilicata','BaData'))
            report.write(getHospitalYear(metadata,'Calabria','ClData'))
            report.write(getHospitalYear(metadata,'Campania','CmData'))
            report.write(getHospitalYear(metadata,'Emilia-Romagna','ErData'))
            report.write(getHospitalYear(metadata,'Friuli-Venezia Giulia','FvData'))
            report.write(getHospitalYear(metadata,'Lazio','LaData'))
            report.write(getHospitalYear(metadata,'Liguria','LiData'))
            report.write(getHospitalYear(metadata,'Lombardia','LoData'))
            report.write(getHospitalYear(metadata,'Marche','MaData'))
            report.write(getHospitalYear(metadata,'Molise','MoData'))
            report.write(getHospitalYear(metadata,'Piemonte','PiData'))
            report.write(getHospitalYear(metadata,'Puglia','PuData'))
            report.write(getHospitalYear(metadata,'Sardegna','SaData'))
            report.write(getHospitalYear(metadata,'Sicilia','SiData'))
            report.write(getHospitalYear(metadata,'Toscana','ToData'))
            report.write(getHospitalYear(metadata,'Trentino-Alto Adige','TaData'))
            report.write(getHospitalYear(metadata,'Umbria','UmData'))
            report.write(getHospitalYear(metadata,"Valle d'Aosta",'VaData'))
            report.write(getHospitalYear(metadata,'Veneto','VeData'))
        else:
            report.write(getRegionSero(metadata))
        report.write(getSTYear(metadata))
        if args.species == "Escherichia coli":
            report.write(getRegionST(metadata))
            report.write(getSero(metadata))
            report.write(getSeroYear(metadata))
            report.write(getStxSub(metadata))
            report.write(getStx12Sero(metadata))
            report.write(getEaeSero(metadata))
        if args.species == "Escherichia coli":
            if metadata.Anno.nunique(dropna = True) == 1:
                insertFile(TOOL_DIR + "/report_tail_stec1.html", report)
            else:
                insertFile(TOOL_DIR + "/report_tail_stec2.html", report)
        else:
            if args.species == "SARS-CoV-2":
                insertFile(TOOL_DIR + "/report_tail_sc2.html", report)
            else:
                insertFile(TOOL_DIR + "/report_tail.html", report)
    finally:
        report.close()

if __name__ == "__main__":
    __main__()


