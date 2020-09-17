#!/bin/bash

wf_ea_script=within_family_ea_simulation.py

########################################
#			WF analyses 
########################################
# attenuated rsquared is ancestry-specific
# EUR ancestry
python3 ${wf_ea_script} \
	--ancestry 'EUR' \
	--main_pheno 'Educational Attainment' \
	--sd_pheno 3.2 \
	--heritability_broad 0.4 \
	--heritability_snp 0.1 \
	--attenuated_rsquared 0.037845 \
	--pleiotropic_trait 'Bipolar Disorder' \
	--heritability_pleio 0.75 \
	--correlation_pleio 0.25 \
	--prevalence_pleio 0.01 &

