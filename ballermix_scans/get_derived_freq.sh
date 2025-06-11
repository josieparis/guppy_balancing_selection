## slurm batch array script

## Load modules
module purge
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

## master directory
MASTER=<dir>

## reference fasta
reference=<dir>

## directory with population-specific VCF files (with AA field)
vcf_files=${MASTER}/<dir>

## output dir
freqs=${MASTER}/outputs/01_raw_derived_freqs

## Set up the chroms
chrs=$(awk '{print $1}' ${reference}.fai | sed "${SLURM_ARRAY_TASK_ID}q;d")

## Set up the pops
#pops=(APHP APLP MHP MLP P GHP GLP TUHP TULP ECHP ECLP)
pops=(MLP P GHP GLP TUHP TULP ECHP ECLP)

for pop in ${pops[@]}
do
vcftools --gzvcf $vcf_files/${pop}.AA.tags.vcf.gz --chr chr${SLURM_ARRAY_TASK_ID} --freq --derived --out ${freqs}/chr${SLURM_ARRAY_TASK_ID}_${pop}.derived
done
