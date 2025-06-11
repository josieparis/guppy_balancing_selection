## Directory containing scripts for the Ballermix analysis

#For B2 DAF first get the derived allele frequencies:

get_derived_freq.sh

#then get genome wide spectrum files for all populations

get_spec.sh

#then run the B2 scan on each population. Here shown just one population as an example

run_ballermix_chrs_APHP_no_win.sh

#For B2 MAF first get the minor allele frequencies:

get_minor_freq.sh

#then get the genome wide spectrum files for all populations

get_spect_maf.sh

#then run the B2 scan for each population. Here shown just one population

run_ballermix_chrs_APHP_no_win_maf.sh

#windowise the Ballermix outputs to compare to other measures

Windowise_ballermix_B2.R
