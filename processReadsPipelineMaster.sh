#!/bin/bash
#PBS -N pipeline
#PBS -l nodes=1:ppn=8
#PBS -o pipeline.out
#PBS -e pipeline.err
#PBS -V
#PBS -q home-yeo
#PBS -W group_list=yeo-group
#PBS -t 1-16
#PBS -l walltime=48:00:00

cmd[1]="/bin/bash processReadsPipeline.sh K562_vitro_N1"
cmd[2]="/bin/bash processReadsPipeline.sh K562_vitro_N2"
cmd[3]="/bin/bash processReadsPipeline.sh K562_vivo_N1"
cmd[4]="/bin/bash processReadsPipeline.sh K562_vivo_N2"

cmd[5]="/bin/bash processReadsPipeline.sh HepG2_vitro_N1"
cmd[6]="/bin/bash processReadsPipeline.sh HepG2_vitro_N2"
cmd[7]="/bin/bash processReadsPipeline.sh HepG2_vivo_N1"
cmd[8]="/bin/bash processReadsPipeline.sh HepG2_vivo_N2"

cmd[9]="/bin/bash processReadsPipeline.sh HeLa_vitro_N1"
cmd[10]="/bin/bash processReadsPipeline.sh HeLa_vitro_N2"
cmd[11]="/bin/bash processReadsPipeline.sh HeLa_vivo_N1"
cmd[12]="/bin/bash processReadsPipeline.sh HeLa_vivo_N2"

cmd[13]="/bin/bash processReadsPipeline.sh vitro_rep1"
cmd[14]="/bin/bash processReadsPipeline.sh vitro_rep2"
cmd[15]="/bin/bash processReadsPipeline.sh vivo_rep1"
cmd[16]="/bin/bash processReadsPipeline.sh vivo_rep2"

eval ${cmd[$PBS_ARRAYID]}
