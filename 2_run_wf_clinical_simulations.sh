#!/bin/bash

wf_script=within_family_simulation.py

################################################################################
#			                WF analyses by ancestry
# condition, heritability, correlation_mz, and prev are same for all ancestries
################################################################################
# EUR ancestry
python3 ${wf_script} \
	--condition 'Type 1 diabetes' 'Type 2 diabetes' 'Breast Cancer (women)' 'Prostate Cancer (men)' \
	'Malignant Melanoma' 'Testicular Cancer (men)' 'Coronary artery disease' 'Hypercholesterolemia' 'Hypertension' \
	'Idiopathic short stature' 'Intellectual disability' \
	--heritability 0.727 0.452 0.31 0.57 0.58 0.37 0.417025 0.44 0.6 0.92 0.55 \
	--correlation_mz 0.727 0.726 0.47 0.57 0.58 0.61 0.417025 0.65 0.67 0.95 0.83 \
 	--prevalence 0.0034 0.3531 0.1286 0.1212 0.0227 0.0041 0.067 0.117 0.46 0.02275 0.02275 \
	--ancestry 'EUR ancestry' \
	--rsquared 0.02940769 0.037922003 0.014808885 0.076605595 0.015842303 0.025845771 0.013730467 0.048492407 0.073421446 0.313 0.052


# AMR ancestry
python3 ${wf_script} \
	--condition 'Type 1 diabetes' 'Type 2 diabetes' 'Breast Cancer (women)' 'Prostate Cancer (men)' \
	'Malignant Melanoma' 'Testicular Cancer (men)' 'Coronary artery disease' 'Hypercholesterolemia' 'Hypertension' \
	'Idiopathic short stature' 'Intellectual disability' \
	--heritability 0.727 0.452 0.31 0.57 0.58 0.37 0.417025 0.44 0.6 0.92 0.55 \
	--correlation_mz 0.727 0.726 0.47 0.57 0.58 0.61 0.417025 0.65 0.67 0.95 0.83 \
 	--prevalence 0.0034 0.3531 0.1286 0.1212 0.0227 0.0041 0.067 0.117 0.46 0.02275 0.02275 \
	--ancestry 'AMR ancestry' \
	--rsquared 0.018379806 0.023701252 0.009255553 0.047878497 0.009901439 0.016153607 0.008581542 0.030307754 0.045888404 0.195625 0.0325


# EAS ancestry
python3 ${wf_script} \
	--condition 'Type 1 diabetes' 'Type 2 diabetes' 'Breast Cancer (women)' 'Prostate Cancer (men)' \
	'Malignant Melanoma' 'Testicular Cancer (men)' 'Coronary artery disease' 'Hypercholesterolemia' 'Hypertension' \
	'Idiopathic short stature' 'Intellectual disability' \
	--heritability 0.727 0.452 0.31 0.57 0.58 0.37 0.417025 0.44 0.6 0.92 0.55 \
	--correlation_mz 0.727 0.726 0.47 0.57 0.58 0.61 0.417025 0.65 0.67 0.95 0.83 \
 	--prevalence 0.0034 0.3531 0.1286 0.1212 0.0227 0.0041 0.067 0.117 0.46 0.02275 0.02275 \
	--ancestry 'EAS ancestry' \
	--rsquared 0.014703845 0.018961002 0.007404443 0.038302797 0.007921151 0.012922886 0.006865233 0.024246203 0.036710723 0.1565 0.026


# AFR ancestry
python3 ${wf_script} \
	--condition 'Type 1 diabetes' 'Type 2 diabetes' 'Breast Cancer (women)' 'Prostate Cancer (men)' \
	'Malignant Melanoma' 'Testicular Cancer (men)' 'Coronary artery disease' 'Hypercholesterolemia' 'Hypertension' \
	'Idiopathic short stature' 'Intellectual disability' \
	--heritability 0.727 0.452 0.31 0.57 0.58 0.37 0.417025 0.44 0.6 0.92 0.55 \
	--correlation_mz 0.727 0.726 0.47 0.57 0.58 0.61 0.417025 0.65 0.67 0.95 0.83 \
 	--prevalence 0.0034 0.3531 0.1286 0.1212 0.0227 0.0041 0.067 0.117 0.46 0.02275 0.02275 \
	--ancestry 'AFR ancestry' \
	--rsquared 0.006535042 0.008427112 0.003290863 0.017023465 0.003520512 0.005743505 0.003051215 0.01077609 0.016315877 0.069555556 0.011555556

