#!/usr/bin/env bash
snakemake -j 100 --cluster-config cluster.json --cluster "sbatch -J test -p {cluster.partition} -N {cluster.node} -cpus-per-task={resources.threads} --mem={resources.mem} -t {resources.time} -o {cluster.stdout} -e {cluster.stderr}"
