#!/bin/bash -l

# # # # # # # # # # # SLURM SETTINGS # # # # # # # # # # # # # #

#SBATCH --account=none # the project account
#SBATCH --job-name=QC_preproc # name of job
#SBATCH --partition=nodes # using CPU nodes
#SBATCH --time=0-04:00:00 # time allocated for job
#SBATCH --mem=4G # memory allocated
#SBATCH --nodes=1 # nodes allocated
#SBATCH --ntasks=1 # number of tasks submitted
#SBATCH --cpus-per-task=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --mail-user=2469736m@student.gla.ac.uk # email address for notifications
#SBATCH --mail-type=END # mail me when my jobs ends
#SBATCH --mail-type=FAIL # mail me if my jobs fails

# # # # # # # # # # # # READ ME # # # # # # # # # # # # # # # # #
: << comment

This script runs fastqc on the files, generates a summary report using multiqc, performs quality trimming with trim galore! and then redoes fastqc and multiqc on the trimmed read files. Also has unused HISAT2 function to remove human reads but i'm not sure if I need to do that yet? 

trim galore settings: 
--phred33 phred encoding, taken from fastqc's prediction (for both 19 and 21!) |  -q 20 phred score cuttoff for quality trimming | -e 0.15 max error rate of 15%, as taken from a consensus of papers | -a forward adaptor sequence | -A reverse adaptor sequence | --max_n 2 maximum number of N bases allowed before the read is discarded. | -o output directory |

cutadapt settings:
-g/-G adaptors are interpreted as being from the 5' of read, or within the 5' region.  and allow  partial matches to 5'. I'm pretty sure the adaptors are from the 5' end, however i did a lot of parameter trial and error using Galaxy and these settings rather than -g/-G ^ADAPTOR, which expects the adaptor only at the start and nowhere within the read yielded the highest quality (note: NOT high, just highest) reads. Also tried polyA tail trimming because there's often this spike of A bases at the 3' end, but that did nothing. These settings at least substantially improve the read quality while not worsening the base bias.
-O minimum overlaps. Default is 3, but since my adaptor searching strategy is more loose-y goose-y i'm increasing it to prevent excessive read loss.
Ran a script before ("alterfilenames.sh") to make the suffixes of I19 and I21 files the same (i.e., remove the "_001" from I21 read file names)


comment
# # # # # # # # # # # # READ ME FIN  # # # # # # # # # # # # # # # 

# activate conda
module load apps/miniforge

# FILEPATHS #

fqpath19=/users/2469736m/msci_project/raw_data/hsct/Ingham19all # path to raw fastq I19 files
fqpath21=/users/2469736m/msci_project/raw_data/hsct/Ingham21all # path to raw fastq I21 files

pretrim_out19='/users/2469736m/msci_project/preprocessing/hsct/Ingham19all/pre_trimQC2' # path to output for I19 pre-trimQC
pretrim_out21='/users/2469736m/msci_project/preprocessing/hsct/Ingham21all/pre_trimQC2' # path to output for I21 pre-trimQC

trimmed_19=/users/2469736m/msci_project/preprocessing/hsct/Ingham19all/post_trim2  # directory to store I19 trimmed reads
trimmed_21=/users/2469736m/msci_project/preprocessing/hsct/Ingham21all/post_trim2  # directory to store I21 trimmed reads

posttrim_out19='/users/2469736m/msci_project/preprocessing/hsct/Ingham19all/post_trim2/post_trimQC2' # directory to store I19 post-trimQC
posttrim_out21='/users/2469736m/msci_project/preprocessing/hsct/Ingham21all/post_trim2/post_trimQC2' # directory to store I21 post-trimQC

# ADAPTORS #

i19_f="CAGCAGCCGCGGTAATAC"
i19_r="CCGTCAATTCCTTTGAGTTT"
# amplifies V4-V5 region

i21_f="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
i21_r="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
# amplifies V3-V4 region amplicon length: ~464bp. Ensure overlap of 480bp so reads can be merged.

# NAMES #
i19='Ingham_2019'
i21='Ingham_2021'
TEST="TEST" # see text area below function calls


# PREPROCESSING FUNCTIONS #

# FASTQC/MULTIQC
fmqc(){
	conda activate MSci
	for sample in $1; do
       		fastqc ${sample} -o $2
	done

	cd $2

	multiqc report 
	multiqc . -n $4
	echo "finished for $5. Output was stored in $2 :)"
	conda deactivate
}

# CUT ADAPT FOR I21 (See readme for explanation)
cut() {
	conda activate trimgalore
	echo "starting quality and adaptor trimming for $1 cutadapt"
	for R1 in ${fqpath21}/*R1.fastq.gz; do
		R2=${R1/R1.fastq.gz/R2.fastq.gz}
		R1o=$(basename ${R1}) # R1 out
		R2o=$(basename ${R2})
		cutadapt -g ${i21_f} -G ${i21_r} -q 25 --max-average-error-rate 0.15 --max-n 1 -O 5 \
			-o $2/${R1o} -p $2/${R2o// /} \
		       	${R1} ${R2} # weird spacing between $3 and R2o which needs to be removed  
	done
	echo "trimming finished for $1. Output was stored in $2 :)"

 	conda deactivate
	echo "starting post-trim quality assessment for $1..."
	fmqc "${trimmed_21}/*.fastq.gz" ${posttrim_out21} 'QC_posttrim21' ${i21}
 
}

# TRIM GALORE! FOR I19 (See readme for explanation)
trimgal() {

	conda activate trimgalore
	echo "starting quality and adaptor trimming for $1 using trim galore!"

        for R1 in ${fqpath19}/*R1.fastq.gz; do
                R2=${R1/R1.fastq.gz/R2.fastq.gz}
                trim_galore --phred33 --paired -a ${i19_f} -a2 ${i19_r} -q 25 -e 0.15 --max_n 1 ${R1} ${R2} \
                        -o $2 
        done

	echo "trimming finished for $1. Ouput was stored in $3 :)"
        conda deactivate

	echo "starting post-trim quality assessment for $1..."
	fmqc "${trimmed_19}/*.fq.gz" ${posttrim_out19} 'QC_posttrim19' ${i19} # in bash, functions can  use variables outwith their scope :D
}



# FUNCTION CALLS #

echo "starting initial quality assessment for ${i19} using fastqc and multiqc..."
#fmqc "${fqpath19}/*" /users/2469736m/msci_project/preprocessing/hsct/Ingham19all/pre_trimQC2  'QC_pretrim19'

echo "initial fastqc and multiqc finished for ${i19}. \n Output was stored in ${pretrim_out19}"
#trimgal ${i19} ${trimmed_19}


echo "starting initial quality assessment for ${i21} using fastqc and multiqc..."
#fmqc "${fqpath21}/*" /users/2469736m/msci_project/preprocessing/hsct/Ingham21all/pre_trimQC2  'QC_pretrim21'

echo "initial fastqc and multiqc finished for ${i21}.\n Ouput was stored in ${pretrim_out21}"
cut ${i21} ${trimmed_21}




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

