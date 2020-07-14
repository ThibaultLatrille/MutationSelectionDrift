#!/usr/bin/env bash
for EXPERIMENT in ./DataEmpirical/Experiments/*Prim*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake build
  snakemake --touch inference
  rm -rf ./CorrelationMatrices
  snakemake --printshellcmds --rerun-incomplete -j 4
  cd ../../..
done

for EXPERIMENT in ./DataEmpirical/PrimatesBinary*; do
  ln -f -s /home/thibault/MutationSelectionDrift/DataEmpirical/run_coevol.sh ${EXPERIMENT}
  ln -f -s /home/thibault/MutationSelectionDrift/DataEmpirical/run_nodetraits.sh ${EXPERIMENT}
  cd ${EXPERIMENT}
  sh ./run_coevol.sh
  sh ./run_nodetraits.sh
  cd ../..
done


for EXPERIMENT in ./DataSimulated/Experiments/df0_*; do
  ln -f -s /home/thibault/MutationSelectionDrift/DataSimulated/run_coevol.sh ${EXPERIMENT}
  cd ${EXPERIMENT}
  # sed -i 's#BURN_IN: 100#BURN_IN: 500#g' config.yaml
  # snakemake --printshellcmds --rerun-incomplete -j 6
  sh ./run_coevol.sh
  cd ../../..
done

for EXPERIMENT in ./DataEmpirical/Experiments/Kappa*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  # sed -i 's#BURN_IN: 300#BURN_IN: 1000#g' config.yaml
  snakemake --unlock
  snakemake build
  snakemake --touch inference
  snakemake --printshellcmds --rerun-incomplete -j 6
  cd ../../..
done

for EXPERIMENT in ./*Isopods*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  sed -i 's#/home/tlatrille/MutationSelectionDrift#/beegfs/data/tlatrill/SimuVsInfer#g' ./*.param ./screen.sh ./snakeslurm.sh
  echo '"#!/usr/bin/env bash"' > snakeslurm.sh
  echo "snakemake -j 99 --cluster 'sbatch -J $(basename $(pwd)) -p long -N 1 -o $(pwd)/slurm.%x.%j.out -e $(pwd)/slurm.%x.%j.err --cpus-per-task={params.threads} --mem={params.mem} -t {params.time}'" >> snakeslurm.sh
  chmod 755 snakeslurm.sh
  sed -i 's#\[node\]#\[node, nodeomega\]#g' config.yaml
  sh ./screen.sh
  cd ../
done

for EXPERIMENT in ./*Mammals*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  sed -i 's#/home/tlatrille/MutationSelectionDrift#/beegfs/data/tlatrill/SimuVsInfer#g' ./*.param ./screen.sh ./snakeslurm.sh
  echo '"#!/usr/bin/env bash"' > snakeslurm.sh
  echo "snakemake -j 99 --cluster 'sbatch -J $(basename $(pwd)) -p long -N 1 -o $(pwd)/slurm.%x.%j.out -e $(pwd)/slurm.%x.%j.err --cpus-per-task={params.threads} --mem={params.mem} -t {params.time}'" >> snakeslurm.sh
  chmod 755 snakeslurm.sh
  sed -i 's#\[node\]#\[node, nodeomega\]#g' config.yaml
  sh ./screen.sh
  cd ../
done

for EXPERIMENT in ./Cat50_OrthoMam_rootedtree*Sample36*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  sed -i 's#RESTART: false#RESTART: true#g' config.yaml
  snakemake --unlock
  sh ./screen.sh
  cd ../
done