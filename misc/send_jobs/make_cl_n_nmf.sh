#! /bin/bash

pipe_dir='misc/send_jobs/pipeline_dir'

pipe_fn="$pipe_dir/pipe_cl.pipe"
rm "$pipe_fn"
touch "$pipe_fn"

# 1 2 3 4 5 6 7 8 9 10 20 30 50 75 100 200 500
for iter in {1..10}
do
  for cl in 150 250 300 350 400 450
  do
  printf "cln_$cl_$iter\t$cl\t-\t.\tn\t2:00:00\t1\t6  \t.\tmodule purge; module load gsl/1.9_64 gsl/2.4 gcc/6.3.0 gmp/6.1.2 R/3.6.0 hdf5/1.10.1; analysis/parameter_benchmarking/cl_n_nmf_job.R $cl $iter\n" >> "$pipe_fn"
  done
done
/home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"
