#!/usr/bin/env
echo "Warming This will only work for python3"

export BLASTDB=/local_uniref/uniref/uniref90


echo "Please input the fasta file of the Beta-barrel protein
      that you want to predict its topology"

read -p "fasta file:" FASTA

DIR="$(dirname "$FASTA")"

#
#
echo "Running PSI-BLAST on $FASTA at $(date)"
time psiblast -query $FASTA -evalue 0.01 -db uniref90.db -num_iterations 3 -out $FASTA.psiblast -out_ascii_pssm  $FASTA.pssm -num_threads 8
echo "Finished running PSI-BLAST on $FASTA at $(date)"
  echo ""
echo "PSI-BLAST done."

echo ""

echo "Running prediction..."

echo ""

cp linear_svmpredictor.py $DIR

python3 linear_svmpredictor.py $FASTA
