#! /bin/bash

pipe_dir='misc/send_jobs/pipeline_dir'

pipe_fn="$pipe_dir/pipe_hvg_nmf.pipe"
rm "$pipe_fn"
touch "$pipe_fn"

# 0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 "uns"
for iter in {1..10}
do
  for hvg in uns
  do
    printf "hvg_$hvg_$iter\t$hvg\t-\t.\tn\t02:00:00\t1\t8\t.\tmodule purge; module load gsl/1.9_64 gsl/2.4 gcc/6.3.0 gmp/6.1.2 R/3.6.0 hdf5/1.10.1; analysis/parameter_benchmarking/hvg_nmf_job.R $hvg $iter\n" >> "$pipe_fn"
  done
done
/home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"

