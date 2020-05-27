#! /bin/bash

pipe_dir='misc/send_jobs/pipeline_dir'

pipe_fn="$pipe_dir/pipe_multiple_tech_nmf.pipe"
rm "$pipe_fn"
touch "$pipe_fn"

for iter in {1..10}
do
  for tech in "Chromium" "inDrop" "C1HT-medium" "C1HT-small" "CEL-Seq2" "ddSEQ" "Drop-Seq" "ICELL8" "MARS-Seq" "Chromium (sn)" "Quartz-Seq2" "SCRB-Seq" "Smart-Seq2"
  do
    printf "$tech_$iter\t$tech\t-\t.\tn\t03:00:00\t1\t4\t.\tmodule purge; module load gsl/1.9_64 gsl/2.4 gcc/6.3.0 gmp/6.1.2 R/3.6.0 hdf5/1.10.1; analysis/parameter_benchmarking/multiple_tech_nmf_job.R $tech $iter\n" >> "$pipe_fn"
  done
done
/home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"
