import argparse as ap
import difflib
from scipy import stats
from collections import Counter

#get samples from index


#polarize to ancestral or major allele freq 
def get_polarized_genotypes(line, index1, index2, index_out):
	
###############
	selected_samples_out = [line[i] for i in index_out]
	samples_genotypes_out = []
	for sample in selected_samples_out:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes_out.append('.')
			samples_genotypes_out.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes_out.append(sample[0])
				samples_genotypes_out.append(sample[0])
			else:
				samples_genotypes_out.append(sample[0])
				samples_genotypes_out.append(sample[1])

###############
	selected_samples1 = [line[i] for i in index1]
	samples_genotypes1 = []
	for sample in selected_samples1:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes1.append('.')
			samples_genotypes1.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes1.append(sample[0])
				samples_genotypes1.append(sample[0])
			else:
				samples_genotypes1.append(sample[0])
				samples_genotypes1.append(sample[1])

###############
	selected_samples2 = [line[i] for i in index2]
	samples_genotypes2 = []
	for sample in selected_samples2:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes2.append('.')
			samples_genotypes2.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes2.append(sample[0])
				samples_genotypes2.append(sample[0])
			else:
				samples_genotypes2.append(sample[0])
				samples_genotypes2.append(sample[1])

	out_geno = [int(i) for i in samples_genotypes_out if i != '.']
	p1_geno = [int(i) for i in samples_genotypes1 if i != '.']
	p2_geno = [int(i) for i in samples_genotypes2 if i != '.']

	if len(set(p1_geno)) == 0:
		ref_allele = int(stats.mode(p2_geno+out_geno)[0]) #not really relevant as our pop will be made of missing data
	elif len(set(out_geno)) == 1 and len(set(p1_geno+p2_geno)) == 2:
		ref_allele = out_geno[0]
	elif len(set(out_geno)) == 1 and len(set(p1_geno+p2_geno)) == 1:
		ref_allele = p1_geno[0]
	elif len(set(out_geno)) == 2 and len(set(p1_geno+p2_geno)) == 1:
		ref_allele = p1_geno[0]
	elif len(set(out_geno)) == 2 and len(set(p1_geno+p2_geno)) == 2:
		ref_allele = int(stats.mode(p1_geno+p2_geno+out_geno)[0]) ### not random! Will always be 0 in case of 50/50
	elif len(set(out_geno)) == 0 and len(set(p1_geno+p2_geno)) > 0:
		ref_allele = int(stats.mode(p1_geno+p2_geno)[0]) ### not random! Will always be 0 in case of 50/50
	elif len(set(out_geno)) == 0 and len(set(p1_geno+p2_geno)) == 0:
		ref_allele = 9
	else:
		ref_allele = 9
		print('WARNING: Non biallelic locus in the vcf. Replaced with missing value for all samples')

	p1_geno_polarized = []
	if ref_allele == 0:
		for i in samples_genotypes1:
			if i == '.':
				p1_geno_polarized.append(9)
			elif int(i) == 0:
				p1_geno_polarized.append(0)
			elif int(i) == 1:
				p1_geno_polarized.append(1)
			elif int(i) == 2:
				p1_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p1_geno_polarized.append(3)###just in case
	elif ref_allele == 1:
		for i in samples_genotypes1:
			if i == '.':
				p1_geno_polarized.append(9)
			elif int(i) == 0:
				p1_geno_polarized.append(1)
			elif int(i) == 1:
				p1_geno_polarized.append(0)
			elif int(i) == 2:
				p1_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p1_geno_polarized.append(3)###just in case
	else:
		for i in samples_genotypes1:
			p1_geno_polarized.append(9)

	return p1_geno_polarized, ref_allele


def replace_all(text):
	table = {'0|0':'1/1','0|1':'1/0','1|0':'0/1','1|1':'0/0'}
	for i, j in table.items():
		text = text.replace(i, j)
	return text


parser = ap.ArgumentParser()
parser.add_argument('-f', '--file', help='An unzipped vcf file', required=True, type=str)
parser.add_argument('-p1', '--popin1', help='Provide txt file with one individual per line', required=True, type=str)
parser.add_argument('-p2', '--popin2', help='Provide txt file with one individual per line', required=True, type=str)
parser.add_argument('-p3', '--popout', help='Provide txt file with one individual per line', required=True, type=str)
args = parser.parse_args()

#define input file
vcf = args.file 

#define output file
outpolarized = open(vcf.replace('.vcf','.pol.vcf'),'w')
# count_polarized_out = 

#get the header
with open(vcf, 'r') as handle:
	for line in handle:
		line = line.rstrip()
		if line.startswith('##'):
			continue
		elif line.startswith('#CHROM'):
			header = line.rstrip().split('\t')
		else:
			break

#get index of ingroup and outgroup samples
samples1 = []
for line in open(args.popin1, 'r').readlines():
    if line != '\n' :
        samples1.append(line.strip())
index_samples1 = []
for sample in samples1:
	index_sample = header.index(sample)
	index_samples1.append(index_sample)

samples2 = []
for line in open(args.popin2, 'r').readlines():
    if line != '\n' :
        samples2.append(line.strip())
index_samples2 = []
for sample in samples2:
	index_sample = header.index(sample)
	index_samples2.append(index_sample)

samples_out = []
for line in open(args.popout, 'r').readlines():
    if line != '\n' :
        samples_out.append(line.strip())

index_samples_out = []
for sample in samples_out:
	index_sample = header.index(sample)
	index_samples_out.append(index_sample)

#go through the vcf and polarize each site
flip_count = []

with open(vcf, 'r') as handle:
	for line in handle:
		if line.startswith('#') or line == '':
			outpolarized.write(line)
		else:
			line = line.rstrip()
			line = line.replace('|','/') # to correct for eventually phased loci as in freebayes MNPs breaking up to SNPs
			line = line.split('\t')
			pos = line[1]
			ref = line[3]
			alt = line[4]
			
			p1_polarized = get_polarized_genotypes(line,index_samples1,index_samples2,index_samples_out)

			anc_allele = p1_polarized[1]
                      #  count_polarized = count(p1_polarized[1])
			if anc_allele == 1:
                                flip = anc_allele,
                                flip_count.append(flip)
                               # print str(len(count))
 				newline = []
				for i in line:
                                        newline.append(replace_all(i.replace('/','|'))) # I have to do this for a weird behaviour of the replace_all function. Don't really understand why it does not work if I try to change '0/0' in '1/1' while it does from '0|0' to '1/1' !!!
                                newline[4] = alt
                                newline[3] = ref
                                outpolarized.write('\t'.join(newline)+'\n')
			elif anc_allele == 0:
                                newline = []
				for i in line:
					newline.append(i) # to correct for the change I made. See above.
                                newline[4] = ref
				newline[4] = ref
				outpolarized.write('\t'.join(newline)+'\n')
			else:
				continue # if there are only missing data ref_allele has been set to 9, the line will be dropped

outpolarized.close()
print "Number of flipped ancestral alleles: ", str(len(flip_count))

