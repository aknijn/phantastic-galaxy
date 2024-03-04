#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import configparser
import os
import mysql.connector
from mysql.connector import errorcode
import pyodbc 

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

class IridaDb:
    def __init__(self, species):
        self._species = species
        if species == 'Escherichia coli' or species == 'Shiga toxin-producing Escherichia coli':
            self._species = 'Escherichia coli'
            self._configFile = TOOL_DIR + '/../phantastic.conf'
            self._allele_header_file = TOOL_DIR + '/data/ecoli.tsv'
        elif species == "Listeria monocytogenes":
            self._configFile = TOOL_DIR + '/../phantastic.conf'
            self._allele_header_file = TOOL_DIR + '/data/listeria.tsv'
        elif species == "Coronavirus" or species == "SARS-CoV-2":
            self._species = 'Coronavirus'
            self._configFile = TOOL_DIR + '/../recovery.conf'
            self._allele_header_file = ''
        else:
            self._configFile = TOOL_DIR + '/../default.conf'
            self._allele_header_file = ''
        config = configparser.ConfigParser()
        config.read(self._configFile)
        dbhost = config['db']['host']
        dbdatabase = config['db']['database']
        dbuser = config['db']['user']
        dbpassword = config['db']['password']
        self._sequence_path = config['fs']['sequence_path']
        self._output_path = config['fs']['output_path']
        config = {
            'user': dbuser, 
            'password': dbpassword, 
            'host': dbhost, 
            'database': dbdatabase
        }
        try:
            self._conn = mysql.connector.connect(**config)
            self._cursor  = self._conn.cursor(buffered=True)
        except mysql.connector.Error as err:
            print(err)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def connection(self):
        return self._conn

    @property
    def cursor(self):
        return self._cursor

    @property
    def species(self):
        return self._species

    @property
    def sequence_path(self):
        return self._sequence_path

    @property
    def output_path(self):
        return self._output_path

    @property
    def allele_header_file(self):
        return self._allele_header_file

    def commit(self):
        self.connection.commit()

    def close(self, commit=True):
        if commit:
            self.commit()
        self.connection.close()

    def execute(self, sql, params=None):
        self.cursor.execute(sql, params or ())

    def fetchall(self):
        return self.cursor.fetchall()

    def fetchone(self):
        return self.cursor.fetchone()

    def query(self, sql, params=None):
        self.cursor.execute(sql, params or ())
        return self.fetchall()

    # PHANtEFSA
    def metadata_for_efsa(self, fileIds):
        if self._species == 'Escherichia coli':
            #sample_id, sample_code, sample_name, fastq1, fastq2, phantastic_type, phantastic_aq, phantastic_seq, phantastic_vir
            sql = "SELECT sample.id as sampleId, metadata_entry.value as sampleCode, sampleName, CONCAT(%s,sf_1.file_path) AS fastq1, \
              CONCAT(%s,sf_2.file_path) AS fastq2, CONCAT(%s,MAX(aof_type.file_path)) AS phantastic_type, CONCAT(%s,MAX(aof_aq.file_path)) \
              AS phantastic_aq,CONCAT(%s,MAX(aof_seq.file_path)) AS phantastic_seq, CONCAT(%s,MAX(aof_vir.file_path)) AS phantastic_vir \
              FROM sample JOIN metadata_entry on(metadata_entry.sample_id = sample.id) JOIN sample_sequencingobject on(sample.id = \
              sample_sequencingobject.sample_id) JOIN sample_sequencingobject AS sso on(sample.id = sso.sample_id) JOIN sequence_file_pair_files \
              AS sfpf on(sso.sequencingobject_id = sfpf.pair_id) JOIN v_sequence_file_pair_files_1 as sfpf_1 on(sfpf.pair_id = sfpf_1.pair_id) \
              JOIN v_sequence_file_pair_files_2 as sfpf_2 on(sfpf.pair_id = sfpf_2.pair_id) JOIN sequence_file AS sf_1 on(sfpf_1.files_id = sf_1.id) \
              JOIN sequence_file AS sf_2 on(sfpf_2.files_id = sf_2.id) JOIN analysis_submission_sequencing_object AS asso on(sfpf.pair_id = \
              asso.sequencing_object_id) JOIN analysis_submission AS asub on(asso.analysis_submission_id = asub.id) JOIN analysis_output_file_map AS \
              aofm_type on(asub.analysis_id = aofm_type.analysis_id) JOIN analysis_output_file AS aof_type on(aofm_type.analysisOutputFilesMap_id = \
              aof_type.id) JOIN analysis_output_file_map AS aofm_aq on(asub.analysis_id = aofm_aq.analysis_id) JOIN analysis_output_file AS aof_aq \
              on(aofm_aq.analysisOutputFilesMap_id = aof_aq.id) JOIN analysis_output_file_map AS aofm_seq on(asub.analysis_id = aofm_seq.analysis_id) \
              JOIN analysis_output_file AS aof_seq on(aofm_seq.analysisOutputFilesMap_id = aof_seq.id) JOIN analysis_output_file_map AS aofm_vir \
              on(asub.analysis_id = aofm_vir.analysis_id) JOIN analysis_output_file AS aof_vir on(aofm_vir.analysisOutputFilesMap_id = aof_vir.id) \
              WHERE metadata_entry.field_id=8 AND aofm_type.analysis_output_file_key = 'phantastic_type' AND aofm_aq.analysis_output_file_key = \
              'phantastic_aq' AND aofm_seq.analysis_output_file_key = 'seqtype' AND aofm_vir.analysis_output_file_key = 'virulotypes' AND \
              sfpf.files_id IN (%s) group by sample.id"
        else:
            sql = "SELECT * FROM sample WHERE sampleName=%s AND sampleName=%s AND sampleName=%s AND sampleName=%s AND sampleName=%s AND sampleName=%s AND id=%s LIMIT 0"
        str_sql = sql % (self.sequence_path, self.sequence_path, self.output_path, self.output_path, self.output_path, self.output_path, fileIds)
        return self.query(str_sql)

    # PHANtEFSA
    def allele_strain(self, sampleCode):
        if self._species == 'Escherichia coli':
            sql = "SELECT allele_strain FROM mlst_ecoli WHERE sample_code=%s"
        elif self._species == "Listeria monocytogenes":
            sql = "SELECT allele_strain FROM mlst_listeria WHERE sample_code=%s"
        else:
            sql = "SELECT * FROM allele_strain WHERE sample_code=%s LIMIT 0"
        return self.query(sql, (sampleCode,))

    # PHANtEFSA
    def update_externalId(self, analyticalPipelineRunId, sampleId):
        if self._species == 'Escherichia coli':
            sql = "UPDATE sample SET externalId=%s WHERE id=%s"
        else:
            sql = "SELECT * FROM sample WHERE externalId=%s AND id=%s LIMIT 0"
        self.execute(sql, (analyticalPipelineRunId, sampleId))
 
    def get_singleend_coverage(self, sample_id):
        sql = "SELECT total_bases/genome_size FROM sample_sequencingobject AS sso \
          INNER JOIN qc_entry  AS so on sso.sequencingObject_id = so.sequencingObject_id \
          INNER JOIN project_sample AS ps on sso.sample_id = ps.sample_id \
          INNER JOIN project AS p on ps.project_id = p.id \
          WHERE genome_size>0 AND sso.sample_id = %s ORDER BY sso.id DESC LIMIT 1"
        self.execute(sql, (sample_id,))
        row = self.fetchone()
        return str(row[0])
        
    # PHANtR
    def metadata_for_report(self, user, fileIds):
        sql_where = " WHERE (email = %s or %s='@iss.it') and files_id in (%s) group by files_id"
        if self._species == "Escherichia coli":
            sql_species = "SELECT ID_ceppo,Anno,Antigen_O,Antigen_H,QC_status,MLST_ST,stx1,stx2,stx_subtype,eae,ehxA,DataCampione,Copertura,CodInterno,vir_tab,amr_tab FROM v_report_ecoli"
            sql = sql_species + sql_where
        elif self._species == "Listeria monocytogenes":
            sql_species = "SELECT ID_ceppo,Anno,QC_status,MLST_ST,MLST_CC,MLST_Lineage,Serogroup,DataCampione,Copertura,CodInterno,vir_tab,amr_tab FROM v_report_listeria"
            sql = sql_species + sql_where
        else:
            sql = "SELECT ID_ceppo,Anno,Antigen_O,Antigen_H,QC_status,MLST_ST,stx1,stx2,stx_subtype,eae,ehxA,DataCampione,Copertura,CodInterno,vir_tab,amr_tab FROM v_report_ecoli " + sql_where + " LIMIT 0"
        str_sql = sql % (user, user[-7:], fileIds)
        return self.query(str_sql)

    # PHANtS
    def metadata_for_summary(self, fileIds):
        if self._species == "Escherichia coli":
            sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_ecoli WHERE files_id IN (%s)")
        elif self._species == "Listeria monocytogenes":
            sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_listeria WHERE files_id IN (%s)")
        elif self._species == "Coronavirus":
            sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_sarscov2 WHERE files_id IN (%s)")
        else:
            sql = ("SELECT Anno,MLST,QC,Regione,Sero,SeroT,Stx12,StxSub,Eae FROM v_summary_ecoli WHERE files_id IN (%s)")
        str_sql = sql % (fileIds)
        return self.query(str_sql)

    # PHANtE
    def metadata_for_export(self, user, fileIds):
        sql_where = " WHERE (lab_users.username = %s or %s='@iss.it') and files_id in (%s) group by files_id"
        if self._species == 'Coronavirus':
            sql_species = "select sample.id AS id,sample.sampleName AS NomeCampione,date_format(sample.collectionDate,'%Y-%m-%d') AS DataCampione,date_format(sample.arrivalDate,'%Y-%m-%d') AS DataArrivoCampione,date_format(sample.createdDate,'%Y-%m-%d') AS DataCaricamento,sample.strain AS CodiceInterno,sample.externalId AS CodiceGISAID,sample.description AS Descrizione,max(case when metadata_entry.field_id = 8 then metadata_entry.value end) AS QC,sample.geographicLocationName AS Regione,sample.geographicLocationName2 AS Provincia,sample.isolationSource AS OrigineIsolato,sample.patientAge AS FasciaEta,sample.patientVaccinationNumber AS NumeroVaccinazioni,sample.patientVaccinationDate AS DataUltimaVaccinazione,sample.collectedBy AS Ospedale,user_group.name AS Laboratorio,f_GetIsolateValue(1,sample.isolate) AS Ospedalizzato,f_GetIsolateValue(2,sample.isolate) AS TerapiaIntensiva,f_GetIsolateValue(3,sample.isolate) AS Reinfetto,f_GetIsolateValue(4,sample.isolate) AS Immunodepresso,f_GetIsolateValue(5,sample.isolate) AS PaeseEsteroAttenzionato,f_GetIsolateValue(6,sample.isolate) AS CampioneCasuale,f_GetIsolateValue(7,sample.isolate) AS CampioneFlash,max(case when metadata_entry.field_id = 7 then metadata_entry.value end) AS Lineage,max(case when metadata_entry.field_id = 54 then metadata_entry.value end) AS Clade,max(case when metadata_entry.field_id = 17 then metadata_entry.value end) AS ORF1ab,max(case when metadata_entry.field_id = 10 then metadata_entry.value end) AS Spike,max(case when metadata_entry.field_id = 9 then metadata_entry.value end) AS ORF3a,max(case when metadata_entry.field_id = 3 then metadata_entry.value end) AS 'E-protein',max(case when metadata_entry.field_id = 16 then metadata_entry.value end) AS 'M-protein',max(case when metadata_entry.field_id = 15 then metadata_entry.value end) AS ORF6,max(case when metadata_entry.field_id = 14 then metadata_entry.value end) AS ORF7a,max(case when metadata_entry.field_id = 13 then metadata_entry.value end) AS ORF7b,max(case when metadata_entry.field_id = 11 then metadata_entry.value end) AS ORF8,max(case when metadata_entry.field_id = 6 then metadata_entry.value end) AS 'N-protein',max(case when metadata_entry.field_id = 18 then metadata_entry.value end) AS ORF10,max(case when metadata_entry.field_id = 35 then metadata_entry.value end) AS Variante,max(case when metadata_entry.field_id = 34 then metadata_entry.value end) AS Sequenza,max(case when metadata_entry.field_id = 52 then if(locate('(',metadata_entry.value) > 0,substr(metadata_entry.value,locate('(',metadata_entry.value) + 1,octet_length(metadata_entry.value) - locate('(',metadata_entry.value) - 2),'') end) AS N_consensus,max(case when qc_entry.total_bases = 0 then 1 else qc_entry.total_bases DIV 29811 end) AS Coverage,lab_users.username AS email,sequence_file_pair_files.files_id AS files_id from metadata_entry join sample on(metadata_entry.sample_id = sample.id) join project_sample on(sample.id = project_sample.sample_id) join user on(user.username = sample.sequencedBy) join user_group_member ugm on(user.id = ugm.user_id) join user_group on(ugm.group_id = user_group.id) join user_group_member lgm on(user_group.id = lgm.group_id) join user lab_users on(lgm.user_id = lab_users.id) join sample_sequencingobject on(sample.id = sample_sequencingobject.sample_id) join sequence_file_pair_files on(sample_sequencingobject.sequencingobject_id = sequence_file_pair_files.pair_id) join qc_entry on(sample_sequencingobject.sequencingobject_id = qc_entry.sequencingObject_id) where project_sample.project_id = 1 and user_group.id > 5 and ugm.role = 'GROUP_MEMBER' and lgm.role = 'GROUP_MEMBER' and "
            sql = sql_species + sql_where
        elif self._species == 'Escherichia coli':
            sql_species = "select sample.id AS id,sample.sampleName AS NomeCampione,date_format(sample.collectionDate,'%Y-%m-%d') AS DataCampione,date_format(sample.arrivalDate,'%Y-%m-%d') AS DataArrivoCampione,date_format(sample.createdDate,'%Y-%m-%d') AS DataCaricamento,sample.strain AS CodiceInterno,sample.externalId AS CodiceEFSA,sample.description AS Descrizione,max(case when metadata_entry.field_id = 13 then metadata_entry.value end) AS QC,sample.geographicLocationName AS Regione,sample.geographicLocationName2 AS Provincia,sample.geographicLocationName3 AS Comune,sample.isolationSource AS OrigineIsolato,sample.collectedBy AS Ospedale,user_group.name AS Laboratorio,sample.isolate AS CondizioneClinica,max(case when metadata_entry.field_id = 3 then metadata_entry.value end) AS Antigen_O,max(case when metadata_entry.field_id = 12 then metadata_entry.value end) AS Antigen_H,max(case when metadata_entry.field_id = 11 then metadata_entry.value end) AS MLST_ST,max(case when metadata_entry.field_id = 14 then metadata_entry.value end) AS stx1,max(case when metadata_entry.field_id = 15 then metadata_entry.value end) AS stx2,max(case when metadata_entry.field_id = 6 then metadata_entry.value end) AS 'stx_subtype',max(case when metadata_entry.field_id = 10 then metadata_entry.value end) AS 'eae',max(case when metadata_entry.field_id = 4 then metadata_entry.value end) AS ehxa,max(case when metadata_entry.field_id = 22 then metadata_entry.value end) AS Virulotipi,max(case when metadata_entry.field_id = 9 then metadata_entry.value end) AS cgMLST_genes_mapped,max(case when metadata_entry.field_id = 7 then metadata_entry.value end) AS Cluster_Id,max(case when qc_entry.total_bases = 0 then 1 else qc_entry.total_bases DIV 5000000 end) AS Coverage,lab_users.username AS email,sequence_file_pair_files.files_id AS files_id from metadata_entry join sample on(metadata_entry.sample_id = sample.id) join user on(user.username = sample.sequencedBy) join user_group_member ugm on(user.id = ugm.user_id) join user_group on(ugm.group_id = user_group.id) join user_group_member lgm on(user_group.id = lgm.group_id) join user lab_users on(lgm.user_id = lab_users.id) join sample_sequencingobject on(sample.id = sample_sequencingobject.sample_id) join sequence_file_pair_files on(sample_sequencingobject.sequencingobject_id = sequence_file_pair_files.pair_id) join qc_entry on(sample_sequencingobject.sequencingobject_id = qc_entry.sequencingObject_id) where ((user_group.id > 5 and ugm.role = 'GROUP_MEMBER' and lgm.role = 'GROUP_MEMBER' and "
            sql = sql_species + sql_where
        elif self._species == 'Listeria monocytogenes':
            sql_species = "select sample.id AS id,sample.sampleName AS NomeCampione,date_format(sample.collectionDate,'%Y-%m-%d') AS DataCampione,date_format(sample.arrivalDate,'%Y-%m-%d') AS DataArrivoCampione,date_format(sample.createdDate,'%Y-%m-%d') AS DataCaricamento,sample.strain AS CodiceInterno,sample.externalId AS CodiceISS,sample.description AS Descrizione,max(case when metadata_entry.field_id = 13 then metadata_entry.value end) AS QC,sample.geographicLocationName AS Regione,sample.geographicLocationName2 AS Provincia,sample.geographicLocationName3 AS Comune,sample.isolationSource AS OrigineIsolato,sample.collectedBy AS Ospedale,user_group.name AS Laboratorio,max(case when metadata_entry.field_id = 19 then metadata_entry.value end) AS Serogroup,max(case when metadata_entry.field_id = 21 then metadata_entry.value end) AS Serotype,max(case when metadata_entry.field_id = 20 then metadata_entry.value end) AS Amplicons,max(case when metadata_entry.field_id = 11 then metadata_entry.value end) AS MLST_ST,max(case when metadata_entry.field_id = 17 then metadata_entry.value end) AS MLST_CC,max(case when metadata_entry.field_id = 18 then metadata_entry.value end) AS 'MLST_Lineage',max(case when metadata_entry.field_id = 9 then metadata_entry.value end) AS 'cgMLST_genes_mapped',max(case when metadata_entry.field_id = 7 then metadata_entry.value end) AS Cluster_Id,max(case when qc_entry.total_bases = 0 then 1 else qc_entry.total_bases DIV 2900000 end) AS Coverage,lab_users.username AS email,sequence_file_pair_files.files_id AS files_id from metadata_entry join sample on(metadata_entry.sample_id = sample.id) join user on(user.username = sample.sequencedBy) join user_group_member ugm on(user.id = ugm.user_id) join user_group on(ugm.group_id = user_group.id) join user_group_member lgm on(user_group.id = lgm.group_id) join user lab_users on(lgm.user_id = lab_users.id) join sample_sequencingobject on(sample.id = sample_sequencingobject.sample_id) join sequence_file_pair_files on(sample_sequencingobject.sequencingobject_id = sequence_file_pair_files.pair_id) join qc_entry on(sample_sequencingobject.sequencingobject_id = qc_entry.sequencingObject_id) where ((user_group.id > 5 and ugm.role = 'GROUP_MEMBER' and lgm.role = 'GROUP_MEMBER' and "
            sql = sql_species + sql_where
        else:
            sql = "select * from announcement WHERE message=%s AND message=%s AND message=%s LIMIT 0"
        str_sql = sql % (user, user[-7:], fileIds)
        return self.query(str_sql)

    # PHANtE
    def header_for_export(self):
        if self._species == 'Coronavirus':
            csv_header = ['id','NomeCampione','DataCampione','DataArrivoCampione','DataCaricamento','CodiceInterno','CodiceGISAID','Descrizione','QC','Regione','Provincia','OrigineIsolato','FasciaEta','NumeroVaccinazioni','DataUltimaVaccinazione','Ospedale','Laboratorio','Ospedalizzato','TerapiaIntensiva','Reinfetto','Immunodepresso','PaeseEsteroAttenzionato','CampioneCasuale','CampioneFlash','Lineage','Clade','ORF1ab','Spike','ORF3a','E-protein','M-protein','ORF6','ORF7a','ORF7b','ORF8','N-protein','ORF10','Variante','Sequenza','N_consensus','Coverage']
        elif self._species == 'Escherichia coli':
            csv_header = ['id','NomeCampione','DataCampione','DataArrivoCampione','DataCaricamento','CodiceInterno','CodiceEFSA','Descrizione','QC','Regione','Provincia','Comune','OrigineIsolato','Ospedale','Laboratorio','CondizioneClinica','Antigen_O','Antigen_H','MLST_ST','stx1','stx2','stx_subtype','eae','ehxa','Virulotipi','cgMLST_genes_mapped','Cluster_Id','Coverage']
        elif self._species == 'Listeria monocytogenes':
            csv_header = ['id','NomeCampione','DataCampione','DataArrivoCampione','DataCaricamento','CodiceInterno','CodiceISS','Descrizione','QC','Regione','Provincia','Comune','OrigineIsolato','Ospedale','Laboratorio','Serogroup','Serotype','Amplicons','MLST_ST','MLST_CC','MLST_Lineage','cgMLST_genes_mapped','Cluster_Id','Coverage']
        else:
            csv_header = []
        return csv_header

    # PHANtEFSA
    def user_in_role(self, username, userrole):
        sql = "SELECT COUNT(*) from user_group_member INNER JOIN user on(user.id=user_group_member.user_id) \
           INNER JOIN user_group on(user_group.id=user_group_member.group_id) WHERE email='%s' and name='%s'"
        str_sql = sql % (username, userrole)
        self.execute(str_sql)
        row = self.fetchone()
        return row[0] > 0

class StecDb:
    def __init__(self, species):
        self._species = species
        if self._species == 'Escherichia coli':
            self._configFile = TOOL_DIR + '/../phantastic.conf'
        elif self._species == "Listeria monocytogenes":
            self._configFile = TOOL_DIR + '/../phantastic.conf'
        elif self._species == "Coronavirus":
            self._configFile = TOOL_DIR + '/../recovery.conf'
        else:
            self._configFile = TOOL_DIR + '/../default.conf'
        config = configparser.ConfigParser()
        config.read(self._configFile)
        dbserver = config['dbmssql']['server']
        dbdatabase = config['dbmssql']['database']
        dbusername = config['dbmssql']['username']
        dbpassword = config['dbmssql']['password']
        self._conn = pyodbc.connect('DRIVER={ODBC Driver 18 for SQL Server};SERVER='+dbserver+';DATABASE='+dbdatabase+';ENCRYPT=yes;UID='+dbusername+';PWD='+ dbpassword)
        self._cursor  = self._conn.cursor(buffered=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def connection(self):
        return self._conn

    @property
    def cursor(self):
        return self._cursor

    @property
    def species(self):
        return self._species

    def commit(self):
        self.connection.commit()

    def close(self, commit=True):
        if commit:
            self.commit()
        self.connection.close()

    def execute(self, sql, params=None):
        self.cursor.execute(sql, params or ())

    def fetchall(self):
        return self.cursor.fetchall()

    def fetchone(self):
        return self.cursor.fetchone()

    def query(self, sql, params=None):
        self.cursor.execute(sql, params or ())
        return self.fetchall()

    def metadata_for_efsa(self, sample_name):
        if self._species == 'Escherichia coli':
            #sample_name, sampId, sampPoint, sampCountry, origCountry, sampArea, sampY, sampM, sampD, sampling_matrix_code, sampling_matrix_free_text, isolId, YEAR(DateOfSampling), MONTH(DateOfSampling), DAY(DateOfSampling)
            sql = "SELECT ISS_ID, sampId, sampPoint, sampCountry, origCountry, sampArea, sampY, sampM, sampD, sampMatCode, '', isolId, \
              YEAR(DateOfSampling), MONTH(DateOfSampling), DAY(DateOfSampling) FROM Samples WHERE ISS_ID=%s"
        else:
            sql = "SELECT * FROM Samples WHERE ISS_ID=%s TOP 0"
        return self.query(sql, (sample_name))
