#! /bin/bash

pipe_dir='misc/send_jobs/pipeline_dir'

pipe_fn="$pipe_dir/pipe_downsampling.pipe"
rm "$pipe_fn"
touch "$pipe_fn"

for iter in {1..10}
do
  for dwnsmpl in 5000 10000 15000 20000 50000
  do
    printf "$dwnsmpl_$iter\t$dwnsmpl\t-\t.\tn\t02:00:00\t1\t6\t.\tmodule purge; module load gsl/1.9_64 gsl/2.4 gcc/6.3.0 gmp/6.1.2 R/3.6.0 hdf5/1.10.1; analysis/parameter_benchmarking/downsample_v3_nmf_job.R $dwnsmpl $iter\n" >> "$pipe_fn"
  done
done
/home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"
