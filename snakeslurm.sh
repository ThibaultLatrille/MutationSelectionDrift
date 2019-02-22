#!/usr/bin/env bash
snakemake -j 100 --cluster-config cluster.json --cluster "sbatch -J SimuInfer -p {cluster.partition} -N {cluster.node} --cpus-per-task={params.threads} --mem={params.mem} -t {params.time} -o {cluster.stdout} -e {cluster.stderr}"
