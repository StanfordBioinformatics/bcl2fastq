#!/bin/sh -e
# set the name of the job
#$ -N bcl2fastq_scg
#
# set the maximum memory usage (per slot)
#$ -l h_vmem=12G
# set the maximum run time
#$ -l h_rt=12:00:00
#
# send mail when job ends or aborts
#$ -m ea
#$ -M pbilling@stanford.edu
#$ -R y
## redirect stoudt
#$ -o /srv/gsfs0/projects/gbsc/workspace/pbilling/code/qsub/logs
#$ -e /srv/gsfs0/projects/gbsc/workspace/pbilling/code/qsub/logs
# check for errors in the job submission options
#$ -w e

# usage: qsub -t 1-8 -A seq_pipeline_cgpm -v R=151117_COOPER_0012_BH3NF7BBXX,B=2,Y=2015,M=dec,C=True run_bcl2fastq_scg.q

module load gbsc/limshostenv/prod
module load /srv/gsfs0/projects/gbsc/workspace/pbilling/modulefiles/scgpm_lims/2.7
module load bcl2fastq2/2.17.1.14

code_dir='/srv/gsfs0/software/gbsc/bcl2fastq'

lims_url='https://uhts.stanford.edu'
lims_token='9af4cc6d83fbfd793fe4'

lane_index=${SGE_TASK_ID}

run_name=${R}
bcl2fastq_version=${B}
year=${Y}
month=${M}
rev_complement=${C} # Default is True, only use False in special cases

if [ ${rev_complement} == 'False' ]
    then
        echo "Not using reverse complement of second barcode or demultiplexing"
        python ${code_dir}/bcl2fastq_scg.py -n -r ${run_name} -l ${lane_index} -b ${bcl2fastq_version} -u ${lims_url} -t ${lims_token} -y ${year} -m ${month}
    elif [ ${rev_complement} == 'True' ] # Default mode: use reverse complement of 2nc barcode
        then 
            python ${code_dir}/bcl2fastq_scg.py -r ${run_name} -l ${lane_index} -b ${bcl2fastq_version} -u ${lims_url} -t ${lims_token} -y ${year} -m ${month}
    else
        echo "Could not determine whether or not to use reverse complement for 2nd barcode. Exiting."
fi
