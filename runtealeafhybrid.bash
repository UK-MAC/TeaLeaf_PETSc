#!/bin/bash

if [ `hostname` == "cy003" ] ; then
  corespernode=24
  decomp="1" # 2 4 6 12 24"
elif [ `hostname` == "cy007" ] ; then
  corespernode=28
  decomp="1 2 4 28"
else
  echo "Please specify cores per node"
  exit 0
fi

benchmarks="tea_bm16_short"
solver_opts="-ksp_type cg -pc_type none"
#solver_opts="-ksp_type cg -pc_type bjacobi -pc_bjacobi_blocks \$((PETSC_NUM_THREADS * SLOTS))"

for nodes in 1 ; do
  for pernode in ${decomp} ; do

    if [ ${pernode} -eq 1 ] || [ ${pernode} -eq ${corespernode} ] ; then
      MPI_DSM="export MPI_DSM_DISTRIBUTE=1"
    else
      MPI_DSM="export MPI_DSM_CPULIST=`seq -s, 0 $((corespernode/pernode)) $((corespernode-1))`"
    fi

    script="runtealeaf_${nodes}_nodes_${pernode}_pernode.pbs"

    cat > ${script} << EOF
#!/bin/bash
#PBS -l walltime=00:40:00
#PBS -N tea_bm16_short
#PBS -j oe
#PBS -l select=${nodes}:ncpus=$((2*corespernode)):mpiprocs=${pernode}:sales_op=none

# Start job from the directory it was submitted
cd \$PBS_O_WORKDIR

echo "Run script:"
echo "---------------------------------------------------------------------"
cat ${script}
echo "---------------------------------------------------------------------"
echo

echo "PBS nodes:"
cat \$PBS_NODEFILE
echo

SLOTS=\$[ \`cat \$PBS_NODEFILE | wc -l\` ]
NODES=\$[ \`uniq \$PBS_NODEFILE | wc -l\` ]
PERNODE=\$((SLOTS / NODES ))

export OMP_NUM_THREADS=\$((${corespernode} / PERNODE))
export PETSC_NUM_THREADS=\$((${corespernode} / PERNODE))

${MPI_DSM}
export MPI_DSM_VERBOSE=1
export KMP_AFFINITY=disabled
export OMPLACE_AFFINITY_COMPAT=ON

# Run benchmarks
for benchmark in ${benchmarks} ; do
  echo
  echo \*\*\* Run benchmark \${benchmark} \*\*\*
  echo
  ln -sf \${benchmark}.in tea.in
  mpiexec_mpt -np \${SLOTS} ./tea_leaf \\
    ${solver_opts} \\
    -threadcomm_type openmp \\
    -threadcomm_nthreads \${PETSC_NUM_THREADS} \\
    -threadcomm_pernode \${PERNODE} \\
    -threadcomm_affinities 0-\$((PETSC_NUM_THREADS * PERNODE)) \\
    -log_summary
  mv tea.out \${benchmark}.out
  echo
  echo \*\*\* Benchmark \${benchmark} completed \*\*\*
  echo
done

echo "Job completed ok!"
EOF

    if [ ! $job ] ; then
      job=`qsub ${script}`
    else
      next_job=`qsub -W depend=afterany:${job} ${script}`
      job=${next_job}
    fi

  done
done
