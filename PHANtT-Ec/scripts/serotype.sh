tooldir="$1";
paired="$2";
fastqfile1="$3";
fastqfile2="$4";
fastafile="$5";

ln -s $fastqfile1 fastqfile1;
ln -s $fastqfile2 fastqfile2;
# FILTER + ASSEMBLE + BLAST FASTQ
if [ $paired = "y" ]
then
  $tooldir/scripts/duk -m filteredO1.fq -k 23 $tooldir/data/O_type.fsa $fastqfile1;
  $tooldir/scripts/duk -m filteredH1.fq -k 23 $tooldir/data/H_type.fsa $fastqfile1;
  cat filteredO1.fq > filteredOH1.fq;
  cat filteredH1.fq >> filteredOH1.fq;
  $tooldir/scripts/duk -m filteredO2.fq -k 23 $tooldir/data/O_type.fsa $fastqfile2;
  $tooldir/scripts/duk -m filteredH2.fq -k 23 $tooldir/data/H_type.fsa $fastqfile2;
  cat filteredO2.fq > filteredOH2.fq;
  cat filteredH2.fq >> filteredOH2.fq;
  $tooldir/scripts/fastq_pair filteredOH1.fq filteredOH2.fq;
  $tooldir/scripts/fastq_pair filteredOH1.fq.single.fq fastqfile2;
  $tooldir/scripts/fastq_pair filteredOH2.fq.single.fq fastqfile1;
  cat filteredOH1.fq.paired.fq > filteredOH1_paired.fq;
  cat filteredOH1.fq.single.fq.paired.fq >> filteredOH1_paired.fq;
  cat fastqfile1.paired.fq >> filteredOH1_paired.fq;
  cat filteredOH2.fq.paired.fq > filteredOH2_paired.fq;
  cat fastqfile2.paired.fq >> filteredOH2_paired.fq;
  cat filteredOH2.fq.single.fq.paired.fq >> filteredOH2_paired.fq;
  dukst1filesize=$(wc -c "filteredOH1_paired.fq" | awk '{print $1}');
  dukst2filesize=$(wc -c "filteredOH2_paired.fq" | awk '{print $1}');
  if [ $dukst1filesize -gt 0 ] && [ $dukst2filesize -gt 0 ]
  then
    perl $tooldir/scripts/spades.pl duk_spades.fasta duk_spades_contig_stats duk_spades_scaffolds duk_spades_scaffold_stats duk_spades_log NODE spades.py --disable-gzip-output --isolate -t \${GALAXY_SLOTS:-16} --pe1-ff --pe1-1 fastq:filteredOH1_paired.fq --pe1-2 fastq:filteredOH2_paired.fq
    rm -r output_dir;
    blastn -query duk_spades.fasta -db $tooldir/data/O_type -task blastn -evalue 0.001 -out duk_O_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;
    blastn -query duk_spades.fasta -db $tooldir/data/H_type -task blastn -evalue 0.001 -out duk_H_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;
  else
    touch duk_O_seqs;
    touch duk_H_seqs;
  fi
else
  $tooldir/scripts/duk -m filteredO1.fq -k 23 $tooldir/data/O_type.fsa $fastqfile1;
  $tooldir/scripts/duk -m filteredH1.fq -k 23 $tooldir/data/H_type.fsa $fastqfile1;
  cat filteredO1.fq > filteredOH1.fq;
  cat filteredH1.fq >> filteredOH1.fq;
  dukstx1filesize=$(wc -c "filteredOH1.fq" | awk '{print $1}');
  if [ $dukstx1filesize -gt 0 ]
  then
    perl $tooldir/scripts/spades.pl duk_spades.fasta duk_spades_contig_stats duk_spades_scaffolds duk_spades_scaffold_stats duk_spades_log NODE spades.py --disable-gzip-output --isolate -t \${GALAXY_SLOTS:-16} --iontorrent -s fastq:filteredOH1.fq;
    rm -r output_dir;
    blastn -query duk_spades.fasta -db $tooldir/data/O_type -task blastn -evalue 0.001 -out duk_O_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;
    blastn -query duk_spades.fasta -db $tooldir/data/H_type -task blastn -evalue 0.001 -out duk_H_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;
  else
    touch duk_O_seqs;
    touch duk_H_seqs;
  fi	
fi
# BLAST FASTA
blastn -query $fastafile -db $tooldir/data/O_type -task blastn -evalue 0.001 -out fasta_O_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;
blastn -query $fastafile -db $tooldir/data/H_type -task blastn -evalue 0.001 -out fasta_H_seqs -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 10 -perc_identity 95.0;

# COMBINE
cat duk_O_seqs > serogroup_O;
cat fasta_O_seqs >> serogroup_O;
cat duk_H_seqs > serogroup_H;
cat fasta_H_seqs >> serogroup_H;
