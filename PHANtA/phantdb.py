import argparse
import configparser
import sys
import os
import mysql.connector
from mysql.connector import errorcode

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

class IridaDb:
    def __init__(self, species):
        self._species = species
        if species == 'Shiga toxin-producing Escherichia coli':
            self._configFile = TOOL_DIR + '/../phantastic.conf'
            self._allele_header_file = TOOL_DIR + '/data/ecoli.tsv'
        elif species == "Listeria monocytogenes":
            self._configFile = TOOL_DIR + '/../phantastic.conf'
            self._allele_header_file = TOOL_DIR + '/data/listeria.tsv'
        elif species == "SARS-CoV-2":
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

    def metadata_for_efsa(self, fileIds):
        if self.species == 'Shiga toxin-producing Escherichia coli':
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

    def allele_strain(self, sampleCode):
        if self.species == 'Shiga toxin-producing Escherichia coli':
            sql = "SELECT allele_strain FROM mlst_ecoli WHERE sample_code=%s"
        elif self.species == "Listeria monocytogenes":
            sql = "SELECT allele_strain FROM mlst_listeria WHERE sample_code=%s"
        else:
            sql = "SELECT * FROM allele_strain WHERE sample_code=%s LIMIT 0"
        return self.query(sql, (sampleCode,))

    def update_externalId(self, analyticalPipelineRunId, sampleId):
        if inspecies == 'Shiga toxin-producing Escherichia coli':
            sql = "UPDATE sample SET externalId=%s WHERE id=%s"
        else:
            sql = "SELECT * FROM sample WHERE externalId=%s AND id=%s LIMIT 0"
        self.execute(sql, (analyticalPipelineRunId, sampleId))
 
    def get_singleend_coverage(self, sample_id):
        sql = "SELECT total_bases/genome_size FROM sequence_file_pair_files AS sfpf \
          INNER JOIN sample_sequencingobject AS sso on sfpf.pair_id = sso.sequencingobject_id \
          INNER JOIN  qc_entry  AS so on sso.sequencingObject_id = so.sequencingObject_id \
          INNER JOIN project_sample AS ps on sso.sample_id = ps.sample_id \
          INNER JOIN project AS p on ps.project_id = p.id \
          WHERE genome_size>0 AND ps.sample_id = %s"
        self.execute(sql, (sample_id,))
        row = self.fetchone()
        return str(row[0])
        
    def user_in_role(self, username, userrole):
        sql = "SELECT COUNT(*) from user_group_member INNER JOIN user on(user.id=user_group_member.user_id) \
           INNER JOIN user_group on(user_group.id=user_group_member.group_id) WHERE email='%s' and name='%s'"
        str_sql = sql % (username, userrole)
        self.execute(str_sql)
        row = self.fetchone()
        return row[0] > 0
