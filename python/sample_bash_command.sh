#! /bin/bash
#$ -cwd
#$ -A xinjunzh
#$ -l h_rt=24:00:00,h_data=8G,h_vmem=16G,arch=intel*
#$ -N misspec-CEUdemo
#$ -o misspec-CEUdemo.output
#$ -j y
#$ -M xinjunzhang@g.ucla.edu
#$ -t 1-100:1

. /u/local/Modules/default/init/modules.sh

module load intel
module load gsl
module load bcftools
module load gcc
module load python

cd /u/scratch/x/xinjunzh/slim_nonAIsweep/misspec/
python3 misspec_CEUdemo.py -s 1 -g $SGE_TASK_ID -d 1
python3 misspec_CEUdemo.py -s 1 -g $SGE_TASK_ID -d 2
python3 misspec_CEUdemo.py -s 1 -g $SGE_TASK_ID -d 3

