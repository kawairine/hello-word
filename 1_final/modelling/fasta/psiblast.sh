export BLASTDB=:/scratch/blastdb

for i in *.fasta
do
  echo "Running PSI-BLAST on $i at $(date)"
  time psiblast -query $i -evalue 0.01 -db uniref90.db -num_iterations 3 -out Output/$i.psiblast -out_ascii_pssm PSSM/$i.pssm -num_threads 8
  echo "Finished running PSI-BLAST on $i at $(date)"
  echo ""
done
echo "Job done."


#export BLASTDB=/local_uniref/uniref/uniref90s
