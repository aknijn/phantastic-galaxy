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
        sql = ("SELECT me_Anno.value AS Anno, me_MLST.value AS MLST, me_QC.value AS QC, me_Reg.value AS Regione, me_AGO.value AS Sero, '-' AS SeroT, me_Cluster.value AS Cluster FROM sequence_file_pair_files " +
              "inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id " +
              "inner join sample_metadata_entry AS sme_Anno on sample_sequencingobject.sample_id = sme_Anno.sample_id " +
              "inner join sample_metadata_entry AS sme_MLST on sample_sequencingobject.sample_id = sme_MLST.sample_id " +
              "inner join sample_metadata_entry AS sme_QC on sample_sequencingobject.sample_id = sme_QC.sample_id " +
              "inner join sample_metadata_entry AS sme_Reg on sample_sequencingobject.sample_id = sme_Reg.sample_id " +
              "inner join sample_metadata_entry AS sme_AGO on sample_sequencingobject.sample_id = sme_AGO.sample_id " +
              "inner join sample_metadata_entry AS sme_Cluster on sample_sequencingobject.sample_id = sme_Cluster.sample_id " +
              "inner join metadata_entry AS me_Anno on sme_Anno.metadata_id = me_Anno.id " +
              "inner join metadata_entry AS me_MLST on sme_MLST.metadata_id = me_MLST.id " +
              "inner join metadata_entry AS me_QC on sme_QC.metadata_id = me_QC.id " +
              "inner join metadata_entry AS me_Reg on sme_Reg.metadata_id = me_Reg.id " +
              "inner join metadata_entry AS me_AGO on sme_AGO.metadata_id = me_AGO.id " +
              "inner join metadata_entry AS me_Cluster on sme_AGO.metadata_id = me_Cluster.id " +
              "WHERE sme_Anno.metadata_KEY = 81 AND sme_MLST.metadata_KEY = 75 AND sme_QC.metadata_KEY = 78 " +
              "AND sme_Reg.metadata_KEY = 69 AND sme_AGO.metadata_KEY = 70 AND sme_Cluster.metadata_KEY = 103 AND sequence_file_pair_files.files_id IN (" + files_id + ")")
    elif species == "Listeria monocytogenes":
        sql = ("SELECT me_Anno.value AS Anno, me_MLST.value AS MLST, me_QC.value AS QC, me_Reg.value AS Regione, me_SG.value AS Sero, me_ST.value AS SeroT, me_Cluster.value AS Cluster FROM sequence_file_pair_files " +
              "inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id " +
              "inner join sample_metadata_entry AS sme_Anno on sample_sequencingobject.sample_id = sme_Anno.sample_id " +
              "inner join sample_metadata_entry AS sme_MLST on sample_sequencingobject.sample_id = sme_MLST.sample_id " +
              "inner join sample_metadata_entry AS sme_QC on sample_sequencingobject.sample_id = sme_QC.sample_id " +
              "inner join sample_metadata_entry AS sme_Reg on sample_sequencingobject.sample_id = sme_Reg.sample_id " +
              "inner join sample_metadata_entry AS sme_SG on sample_sequencingobject.sample_id = sme_SG.sample_id " +
              "inner join sample_metadata_entry AS sme_ST on sample_sequencingobject.sample_id = sme_ST.sample_id " +
              "inner join sample_metadata_entry AS sme_Cluster on sample_sequencingobject.sample_id = sme_Cluster.sample_id " +
              "inner join metadata_entry AS me_Anno on sme_Anno.metadata_id = me_Anno.id " +
              "inner join metadata_entry AS me_MLST on sme_MLST.metadata_id = me_MLST.id " +
              "inner join metadata_entry AS me_QC on sme_QC.metadata_id = me_QC.id " +
              "inner join metadata_entry AS me_Reg on sme_Reg.metadata_id = me_Reg.id " +
              "inner join metadata_entry AS me_SG on sme_SG.metadata_id = me_SG.id " +
              "inner join metadata_entry AS me_ST on sme_ST.metadata_id = me_ST.id " +
              "inner join metadata_entry AS me_Cluster on sme_Cluster.metadata_id = me_Cluster.id " +
              "WHERE sme_Anno.metadata_KEY = 81 AND sme_MLST.metadata_KEY = 75 AND sme_QC.metadata_KEY = 78 " +
              "AND sme_Reg.metadata_KEY = 69 AND sme_SG.metadata_KEY = 63 AND sme_ST.metadata_KEY = 96 AND sme_Cluster.metadata_KEY = 103 AND sequence_file_pair_files.files_id IN (" + files_id + ")")
    elif species == "SARS-CoV-2":
        sql = ("SELECT me_Anno.value AS Anno, LEFT(me_MLST.value,INSTR(me_MLST.value, '(')-2) AS MLST, me_QC.value AS QC, me_Reg.value AS Regione, '-' AS Sero, '-' AS SeroT, '-' AS Cluster FROM sequence_file_pair_files " +
              "inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id " +
              "inner join sample_metadata_entry AS sme_Anno on sample_sequencingobject.sample_id = sme_Anno.sample_id " +
              "inner join sample_metadata_entry AS sme_MLST on sample_sequencingobject.sample_id = sme_MLST.sample_id " +
              "inner join sample_metadata_entry AS sme_QC on sample_sequencingobject.sample_id = sme_QC.sample_id " +
              "inner join sample_metadata_entry AS sme_Reg on sample_sequencingobject.sample_id = sme_Reg.sample_id " +
              "inner join metadata_entry AS me_Anno on sme_Anno.metadata_id = me_Anno.id " +
              "inner join metadata_entry AS me_MLST on sme_MLST.metadata_id = me_MLST.id " +
              "inner join metadata_entry AS me_QC on sme_QC.metadata_id = me_QC.id " +
              "inner join metadata_entry AS me_Reg on sme_Reg.metadata_id = me_Reg.id " +
              "WHERE sme_Anno.metadata_KEY = 81 AND sme_MLST.metadata_KEY = 142 AND sme_QC.metadata_KEY = 78 " +
              "AND sme_Reg.metadata_KEY = 69 AND sequence_file_pair_files.files_id IN (" + files_id + ")")
    else:
        sql = ("SELECT me_Anno.value AS Anno, me_MLST.value AS MLST, me_QC.value AS QC, me_Reg.value AS Regione, '-' AS Sero, '-' AS SeroT, me_Cluster.value AS Cluster FROM sequence_file_pair_files " +
              "inner join sample_sequencingobject on sequence_file_pair_files.pair_id = sample_sequencingobject.sequencingobject_id " +
              "inner join sample_metadata_entry AS sme_Anno on sample_sequencingobject.sample_id = sme_Anno.sample_id " +
              "inner join sample_metadata_entry AS sme_MLST on sample_sequencingobject.sample_id = sme_MLST.sample_id " +
              "inner join sample_metadata_entry AS sme_QC on sample_sequencingobject.sample_id = sme_QC.sample_id " +
              "inner join sample_metadata_entry AS sme_Reg on sample_sequencingobject.sample_id = sme_Reg.sample_id " +
              "inner join sample_metadata_entry AS sme_SG on sample_sequencingobject.sample_id = sme_SG.sample_id " +
              "inner join sample_metadata_entry AS sme_Cluster on sample_sequencingobject.sample_id = sme_Cluster.sample_id " +
              "inner join metadata_entry AS me_Anno on sme_Anno.metadata_id = me_Anno.id " +
              "inner join metadata_entry AS me_MLST on sme_MLST.metadata_id = me_MLST.id " +
              "inner join metadata_entry AS me_QC on sme_QC.metadata_id = me_QC.id " +
              "inner join metadata_entry AS me_Reg on sme_Reg.metadata_id = me_Reg.id " +
              "inner join metadata_entry AS me_Cluster on sme_AGO.metadata_id = me_Cluster.id " +
              "WHERE sme_Anno.metadata_KEY = 81 AND sme_MLST.metadata_KEY = 75 AND sme_QC.metadata_KEY = 78 " +
              "AND sme_Reg.metadata_KEY = 69 AND sme_Cluster.metadata_KEY = 103 AND sequence_file_pair_files.files_id IN (" + files_id + ")")

    try:
        cnx = mysql.connector.connect(**config)
        cursor = cnx.cursor(buffered=True)
        cursor.execute(sql)
        records = pd.DataFrame(cursor.fetchall(),columns=['Anno','MLST','QC','Regione','Sero','SeroT','Cluster'])
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

def getYearRegion(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Regione"], columns='Anno', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    dfregions = dfpivot.index
    strYR = "var AnnoRegioneData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lstregions = []
    for idx, row in dfpivot.iterrows():
        region = idx.replace("__sq__", "'")
        lstregiondata = []
        for dfkey in dfkeys:
            lstregiondata.append(row[dfkey])
        strregiondata = ",".join([str(item) for item in lstregiondata])
        lstregions.append("{label: \"" + region + "\", backgroundColor : getRCH(), data: [" + strregiondata + "]}")
    strYR = strYR + ",".join(lstregions) + "]};\n"
    return strYR

def getST(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], values='QC', aggfunc=len, fill_value=0)
    dfSTs = dfpivot.index
    strSR = "var STData = {labels: [" + ",".join(["'" + item + "'" for item in dfSTs]) + "],datasets: [{backgroundColor: ["
    strSR = strSR + ",".join(["getRCH()" for item in dfSTs]) + "], data: ["
    lstSTdata = []
    for idx, row in dfpivot.iterrows():
        lstSTdata.append(str(row[0]))
    strSR = strSR + ",".join(lstSTdata) + "]}]};\n"
    return strSR

def getRegionSero(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["Sero"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    dfseros = dfpivot.index
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

def getRegionLineage(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], columns='Regione', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    dfseros = dfpivot.index
    strRL = "var RegioneLineageData = {labels: [" + ",".join(["\"" + item.replace("__sq__", "'") + "\"" for item in dfkeys]) + "],datasets: ["
    lstlineages = []
    for idx, row in dfpivot.iterrows():
        sero = idx
        lstlineagedata = []
        for dfkey in dfkeys:
            lstlineagedata.append(row[dfkey])
        strlineagedata = ",".join([str(item) for item in lstlineagedata])
        lstlineages.append("{label: \"" + sero + "\", backgroundColor : getRCH(), data: [" + strlineagedata + "]}")
    strRL = strRL + ",".join(lstlineages) + "]};\n"
    return strRL

def getYearST(dataframe):
    dfpivot = pd.pivot_table(dataframe,index=["MLST"], columns='Anno', values='QC', aggfunc=len, fill_value=0)
    dfkeys = dfpivot.keys().tolist()
    dfsts = dfpivot.index
    strYS = "var AnnoSTData = {labels: [" + ",".join(["'" + item + "'" for item in dfkeys]) + "],datasets: ["
    lststs = []
    for idx, row in dfpivot.iterrows():
        lststdata = []
        for dfkey in dfkeys:
            lststdata.append(row[dfkey])
        strstdata = ",".join([str(item) for item in lststdata])
        lststs.append("{label: \"" + idx + "\", backgroundColor : getRCH(), data: [" + strstdata + "]}")
    strYS = strYS + ",".join(lststs) + "]};\n"
    return strYS

def getClusterTable(dataframe):
    df_clusters = dataframe.groupby('Cluster')
    #Column Count
    df_count = df_clusters.count()[['MLST']]
    df_count.columns = ['Numero di casi']
    #Columns Sequence type & Serogroup
    df_first = df_clusters.first()
    df_first.drop(['QC'],axis=1,inplace=True)
    df_first.drop(['Anno'],axis=1,inplace=True)
    df_first.drop(['Regione'],axis=1,inplace=True)
    df_first = df_first[['MLST','Sero','SeroT']]
    df_first.columns = ['Sequence type', 'Serogruppo', 'Serotipo']
    #Column Years
    df_years = dataframe.copy()
    df_years.drop(['QC'],axis=1,inplace=True)
    df_years.drop(['MLST'],axis=1,inplace=True)
    df_years.drop(['Regione'],axis=1,inplace=True)
    df_years.drop(['Sero'],axis=1,inplace=True)
    df_years.drop(['SeroT'],axis=1,inplace=True)
    df_years = df_years.drop_duplicates()
    df_years = df_years.sort_values(by=['Anno']).groupby('Cluster')
    df_years = df_years['Anno'].apply(lambda tags: ', '.join(tags))
    #Column Regions
    df_regions = dataframe.copy()
    df_regions.drop(['QC'],axis=1,inplace=True)
    df_regions.drop(['MLST'],axis=1,inplace=True)
    df_regions.drop(['Anno'],axis=1,inplace=True)
    df_regions.drop(['Sero'],axis=1,inplace=True)
    df_regions.drop(['SeroT'],axis=1,inplace=True)
    df_regions = df_regions.drop_duplicates()
    df_regions = df_regions.sort_values(by=['Regione']).groupby('Cluster')
    df_regions = df_regions['Regione'].apply(lambda tags: ', '.join(tags))
    #df_regions.columns = ['Cluster', 'Regioni coinvolti']
    #Merge dataframes
    df_intermediate1 = df_count.merge(df_first, left_on='Cluster', right_on='Cluster')
    df_intermediate2 = df_intermediate1.merge(df_years, left_on='Cluster', right_on='Cluster')
    df_final = df_intermediate2.merge(df_regions, left_on='Cluster', right_on='Cluster')
    df_final = df_final.filter(regex ='^[^_]*_[^_]*$', axis=0)
    strCT = df_final.to_html(classes='table table-cross', border=0)
    return strCT  

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='phants filenames')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--phants_stat', dest='phants_stat', help='phants stat report html file')
    args = parser.parse_args()
    if args.species == "Shiga toxin-producing Escherichia coli":
        args.species = "Escherichia coli"
    if args.species == "Coronavirus":
        args.species = "SARS-CoV-2"
    metadata = getMetadata(args.input_files, args.species)
    try:
        report = open(args.phants_stat, 'w')
        # write head html
        if args.species == "SARS-CoV-2":
            insertFile(TOOL_DIR + "/report_head_sc2.html", report)
        else:
            insertFile(TOOL_DIR + "/report_head.html", report)
        if args.species != "SARS-CoV-2":
            report.write(getClusterTable(metadata))
        report.write("</td></tr></table>\n")
        report.write("<div>ISS: Riepilogo per <i>%s</i>, elaborato %s</div>\n<script>\n" % (args.species, datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")))
        report.write(getPassedFailed(metadata))
        report.write(getYearRegion(metadata))
        report.write(getST(metadata))
        if args.species == "SARS-CoV-2":
            report.write(getRegionLineage(metadata))
        else:
            report.write(getRegionSero(metadata))
        report.write(getYearST(metadata))
        if args.species == "SARS-CoV-2":
            insertFile(TOOL_DIR + "/report_tail_sc2.html", report)
        else:
            insertFile(TOOL_DIR + "/report_tail.html", report)
    finally:
        report.close()

if __name__ == "__main__":
    __main__()


