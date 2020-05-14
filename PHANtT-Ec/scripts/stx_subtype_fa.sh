tooldir="$1";
fastafile="$2";

# SEQUENCE SEARCH
blastn -query $fastafile -db $tooldir/data/stx -task blastn -evalue 0.001 -out shigatoxin -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 1 -perc_identity 95.0;
# SHIGATOXINTYPER: FILTER, CUT AND CONCATENATE SEQUENCE SEARCH OUTPUT
awk -F '\t' '($3>95 && $4>1200) { print $2 FS $3 FS $4 FS $16 }' shigatoxin > shigatoxin_fc;