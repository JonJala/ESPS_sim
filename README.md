# Embryo Selection simulation 
This repo includes the code used for between- and within-family simulations for educational attainment and clinical phenotypes, as described in Challenges With Embryo Selection Using Polygenic Scores (Turley et al.). 

The pieces of code are: 
* `between_family_ea_simulation.py` -- This is used to for the between-family EA simulation. It can be run with `1_run_bf_ea_simulation.sh`. **This code produces data used for Figure 1.**
* `within_family_ea_simulation.py` -- This is used for the within-family EA simulation. It can be run with `1_run_wf_ea_simulation.sh`. **This code produces data used for Figure 1 and values for the pleiotropy analysis of EA and bipolar disorder.**
* `within_family_simulation.py` -- This is used for the within-family clinical phenotypes simulation. Simulations for type 1 diabetes, type 2 diabetes, breast/prostate/testicular cancers, malignant melanoma, coronary artery disease, hypercholesterolemia, hypertension, idiopathic short stature, and intellectual disability can be run using `2_run_wf_clinical_simulations.sh`. Note that clinical phenotypes must be binary (e.g. hypertension is binary, whereas blood pressure is continuous). **This code produces data that is used to create Figure 2 and all supplementary figures and tables.**
