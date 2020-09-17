"""
Purpose: To simulate expected educational attainment gains from embryo selection between families.
Date: 10/09/2019
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from between_family_ea_simulation import (
	get_random_index, 
	get_max_pgs_index,
	select_embryos_by_index,
	calc_phenotype_diffs
	)
from scipy.stats import norm
import argparse



def calc_within_family_values(n, num_embryos, heritability, correlation_mz, rsquared):
    """
    Purpose: To get the ghat_i and y_i for each family pair, where i={1,...,num_embryos}.

    Arguments:
        n: integer number of parent pairs
        num_embryos: integer number of embryos for each parent pair
        heritability: heritability of clinical trait
        correlation_mz: twin correlation of clinical trait
        rsquared: Ancestry-specific R^2 value for PGS prediction of trait

    Returns:
        {'pgs':df_pgs, 'liability':df_liability}. 
        Each dataframe has size (n x num_embryos) and holds polygenic scores and phenotype liability values, respectively.
    """
    df_a = pd.DataFrame()
    df_pgs = pd.DataFrame()
    df_liability = pd.DataFrame()

    # calculate average parent additive component
    a_mean = np.random.normal(loc=0, scale=(heritability/2)**0.5, size=int(n))

    # calculate shared environmental component (constant within-family)
    c = np.random.normal(loc=0, scale=(correlation_mz - heritability)**0.5, size=int(n))

    # calculate liability and PGS for clinical trait, for each embryo
    for i in range(num_embryos):
    	# generate individual embryo additive components
    	df_a['a_{}'.format(i)] = np.random.normal(loc=a_mean, scale=(heritability/2)**0.5, size=int(n))
        # calculate PGS
    	df_pgs['ghat_{}'.format(i)] = np.random.normal(loc=df_a['a_{}'.format(i)], \
            scale=(heritability**2/rsquared - heritability)**0.5, size=int(n))
        # calculate phenotypic liability
    	df_liability['liability_{}'.format(i)] =  np.random.normal(loc=(df_a['a_{}'.format(i)] + c), \
    		scale=(1 - correlation_mz)**0.5, size=int(n))

    return {'pgs':df_pgs, 'liability':df_liability}


def process_arguments():
    """
    Parses command line arguments.

    Args:
    -----
    None

    Returns:
    --------
    parser: :class: argparse.Namespace
        arguments passed in from command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=1000000, type=int, 
                        help='Number of parent pairs to simulate. Defaults to 1e6.')
    parser.add_argument('--embryos', default=10, type=int, 
                        help='Number of embryos from which to choose. Defaults to 10.')
    parser.add_argument('--ancestry', required=True, type=str, 
    					help='Ancestry of interest.')
    parser.add_argument('--heritability', required=True, type=float, nargs='+',
    					help='List of heritabilities on the liability scale for conditions. \
    					Index of heritabilities must match index of rsquared, correlation_mz, and prevalence.')
    parser.add_argument('--rsquared', required=True, type=float, nargs='+',
    					help='List of r2 for conditions. Index of r2 must match index of \
    					heritability, correlation_mz, and prevalence.')
    parser.add_argument('--correlation_mz', required=True, type=float, nargs='+',
    					help='List of monozygotic twin correlations on the liability scale for conditions. \
    					Index of correlations must match index of rsquared, heritabilities, and prevalence.')
    parser.add_argument('--prevalence', required=True, type=float, nargs='+',
    					help='List of prevalences for conditions. \
    					Index of prevalences must match index of rsquared, heritabilities, and correlations.')
    parser.add_argument('--condition', required=True, type=str, nargs='+',
    					help='Name of conditions. \
    					Index must match index of rsquared, heritabilities, and correlations.')
    return parser.parse_args()


if __name__ == "__main__":
	# import arguments
	args = process_arguments()
	N = args.n 
	NUM_EMBRYOS = args.embryos
	HERITABILITY = args.heritability # list
	RSQUARED = args.rsquared # list
	CORR_MZ = args.correlation_mz # list
	PREVALENCE = args.prevalence # list
	CONDITION = args.condition # list

	# check lengths
	assert len(HERITABILITY) == len(RSQUARED), 'Your lists aren\'t the same length!'
	assert len(HERITABILITY) == len(CORR_MZ), 'Your lists aren\'t the same length!'
	assert len(HERITABILITY) == len(PREVALENCE), 'Your lists aren\'t the same length!'
	assert len(HERITABILITY) == len(CONDITION), 'Your lists aren\'t the same length!'

	# print intro
	print('This is a simulation for within-family selection of embryos.')
	print('This analysis is for parents of ' + args.ancestry + ' ancestry who are choosing from ' + str(NUM_EMBRYOS) + ' embryos.')

	#### begin simualation
	for i in range(len(args.condition)):
		# calculate values
		values = calc_within_family_values(n=N, num_embryos=NUM_EMBRYOS, 
			heritability=HERITABILITY[i], correlation_mz=CORR_MZ[i], rsquared=RSQUARED[i])

		# generate indices
		rand_index = get_random_index(NUM_EMBRYOS, N)
		max_index = get_max_pgs_index(values['pgs'])

		# get max/random liability values
		max_liability = select_embryos_by_index(values['liability'], max_index)
		rand_liability = select_embryos_by_index(values['liability'], rand_index)

		# get fraction of individuals who will have disease (convert liability to binary trait)
		max_frac = (max_liability <= norm.ppf(PREVALENCE[i])).astype(int).mean()
		rand_frac = (rand_liability <= norm.ppf(PREVALENCE[i])).astype(int).mean()

		# print summary
		print('For ' + CONDITION[i] + ', the within-family prevalence is ' + str(rand_frac['selected_values']) + \
				' in random embryos and ' + str(max_frac['selected_values']) + ' in selected embryos.')

