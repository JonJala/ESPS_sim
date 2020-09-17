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
import argparse


def calc_within_family_values(n, no_embryos, hsquared_broad, hsquared_snp, rsquared):
    """
    Purpose: To get the PGS and unscaled EA phenotype value for all embryos for each parent pair. 

    Arguments:
        n: integer number of parent pairs
        no_embryos: integer number of embryos for each parent pair
        hsquared_broad: broad heritability
        hsquared_snp: SNP heritability
        rsquared: Ancestry-specific R^2 value for PGS prediction of trait

    Returns:
        {'ghat':df_ghat, 'y':df_a}. 
        Each dataframe has size (n x no_embryos) and holds polygenic scores and phenotype values (unscaled), respectively.
    """
    df_g = pd.DataFrame()
    df_xi = pd.DataFrame()
    df_a = pd.DataFrame()
    df_ghat = pd.DataFrame()
    df_pleio = pd.DataFrame()

    # calc average parent g for each parent pair
    g_mean = np.random.normal(loc=0, scale=(hsquared_snp/2)**0.5, size=int(n))

    # calc average parent xi for each parent pair
    xi_mean = np.random.normal(loc=0, scale=((hsquared_broad-hsquared_snp)/2)**0.5, size=int(n))

    # generate EA PGS and EA phenotyp values for each embryo. Note: g and xi are used to calculate a and ghat.
    for i in range(no_embryos):
        df_g['g_{}'.format(i)] = np.random.normal(loc=g_mean, scale=(hsquared_snp/2)**0.5, size=int(n))
        df_xi['xi_{}'.format(i)] = np.random.normal(loc=xi_mean, \
            scale=((hsquared_broad-hsquared_snp)/2)**0.5, size=int(n))
        df_a['a_{}'.format(i)] = df_g['g_{}'.format(i)] + df_xi['xi_{}'.format(i)]
        df_ghat['ghat_{}'.format(i)] = np.random.normal(loc=df_g['g_{}'.format(i)], \
            scale=(hsquared_snp**2/rsquared - hsquared_snp)**0.5, size=int(n))

    return {'ghat':df_ghat, 'y':df_a}


def calc_within_family_values_and_pleiotropy(n, no_embryos, hsquared_broad, hsquared_snp, rsquared, r_g, hsquared_pleio):
    """
    Purpose: To get the PGS and unscaled EA phenotype value for all embryos for each parent pair. 
    This function also returns the phenotype values for pleiotropic trait. 

    Arguments:
        n: integer number of parent pairs
        no_embryos: integer number of embryos for each parent pair
        hsquared_broad: broad heritability
        hsquared_snp: SNP heritability
        rsquared: Ancestry-specific R^2 value for PGS prediction of trait
        r_g: correlation coefficient between EA and the pleiotropic trait
        hsquared_pleio: heritability of pleiotropic trait

    Returns:
        {'ghat':df_ghat, 'y':df_a, 'pleio':df_pleio}. 
        Each dataframe has size (n x no_embryos) and holds EA polygenic scores, EA phenotype values (unscaled), and pleiotropic phenotype values respectively.
    """
    df_g = pd.DataFrame()
    df_xi = pd.DataFrame()
    df_a = pd.DataFrame()
    df_ghat = pd.DataFrame()
    df_pleio = pd.DataFrame()

    # calc average parent g for each parent pair
    g_mean = np.random.normal(loc=0, scale=(hsquared_snp/2)**0.5, size=int(n))

    # calc average parent xi for each parent pair
    xi_mean = np.random.normal(loc=0, scale=((hsquared_broad-hsquared_snp)/2)**0.5, size=int(n))

    # generate EA PGS, EA phenotyp values, and pleiotropic phenotype values for each embryo. Note: g and xi are used to calculate a and ghat.
    for i in range(no_embryos):
        df_g['g_{}'.format(i)] = np.random.normal(loc=g_mean, scale=(hsquared_snp/2)**0.5, size=int(n))
        df_xi['xi_{}'.format(i)] = np.random.normal(loc=xi_mean, \
            scale=((hsquared_broad-hsquared_snp)/2)**0.5, size=int(n))
        df_a['a_{}'.format(i)] = df_g['g_{}'.format(i)] + df_xi['xi_{}'.format(i)]
        df_ghat['ghat_{}'.format(i)] = np.random.normal(loc=df_g['g_{}'.format(i)], \
            scale=(hsquared_snp**2/rsquared - hsquared_snp)**0.5, size=int(n))
        df_pleio['a_{}'.format(i)] = np.random.normal(loc=r_g*((hsquared_pleio/hsquared_broad)**0.5)*df_a['a_{}'.format(i)], \
            scale=(hsquared_pleio*(1-r_g**2))**0.5, size=int(n))

    return {'ghat':df_ghat, 'y':df_a, 'pleio':df_pleio}


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
    parser.add_argument('--sd_pheno', required=True, type=float, 
                        help='Standard deviation of the main phenotype of interest.')
    parser.add_argument('--heritability_broad', required=True, type=float, 
                        help='Broad sense heritability for the main phenotype of interest.')
    parser.add_argument('--heritability_snp', required=True, type=float,
                        help='SNP-based heritability for the main phenotype of interest. Must be less than or equal to heritability_broad.')
    parser.add_argument('--attenuated_rsquared', required=True, type=float,
                        help='Attenuated R^2 for main phenotype of interest.')
    parser.add_argument('--heritability_pleio', default=None, type=float, 
                        help='Heritability for pleiotropic trait. Defaults to None.')
    parser.add_argument('--correlation_pleio', default=None, type=float, 
                        help='Correlation between main trait and pleiotropic trait. Defaults to None.')
    parser.add_argument('--prevalence_pleio', default=None, type=float,
                        help='Prevalence of pleiotropic in the general population. Defaults to None.')
    parser.add_argument('--pleiotropic_trait', default=None, type=str,
                        help='Name of pleiotropic trait. Defaults to None. Trait must be binary, not continuous.')
    parser.add_argument('--main_pheno', default='Educational Attainment', type=str, 
                        help='Name of main trait of interest. Main trait must be continuous. Defaults to "Educational Attainment".')
    return parser.parse_args()


if __name__ == "__main__":
    # import arguments
    args = process_arguments()

    # check broad heritability >= snp heritability
    assert args.heritability_snp <= args.heritability_broad, 'SNP heritability must be less than or equal to broad-sense heritability!'

    # if no pleiotropic trait is specified, then do this version
    if args.heritability_pleio==None or args.correlation_pleio==None or args.prevalence_pleio==None or args.pleiotropic_trait==None:
        # print intro
        print('This is a within-family simulation for embryo selection using the polygenic score for ' + args.main_pheno + '.')
        print('This analysis is done for ' + str(args.n) + ' parents of ' + args.ancestry + ' ancestry choosing between ' + str(args.embryos) + ' embryos.')
        print('Note: there is no pleiotropic trait considered in this simulation.')

        # calculate within-family values
        values = calc_within_family_values(n=args.n, no_embryos=args.embryos, hsquared_broad=args.heritability_broad, 
            hsquared_snp=args.heritability_snp, rsquared=args.attenuated_rsquared)

        # generate indices
        rand_index = get_random_index(args.embryos, args.n)
        max_index = get_max_pgs_index(values['ghat'])

        # get max y values
        max_y = select_embryos_by_index(values['y'], max_index)

        # get random y values
        rand_y = select_embryos_by_index(values['y'], rand_index)

        # calculate difference
        diffs_y = calc_phenotype_diffs(max_y, rand_y, args.sd_pheno)
        mean = diffs_y['selected_values'].mean()
        interval = diffs_y['selected_values'].std() * 1.96

        # print results
        print('For parents of ' + str(args.ancestry) + ' ancestry, the within-family 95-percent prediction interval for the difference in ' + args.main_pheno +  \
            ' between random and selected embryos is ' + '%.2f' % mean + ' +/- ' +'%.2f' % interval + '.')
        print('End of analysis.')

    # if pleitropic trait is specified, then do this version
    else:
        # print intro
        print('This is a within-family simulation for embryo selection using the polygenic score for ' + args.main_pheno + '.')
        print('This analysis is done for ' + str(args.n) + ' parents of ' + args.ancestry + ' ancestry choosing between ' + str(args.embryos) + ' embryos.')
        print('Note: we are considering the pleiotropic trait of ' + args.pleiotropic_trait + ' in this analysis.')
        
        # calculate within-family values
        values = calc_within_family_values_and_pleiotropy(n=args.n, no_embryos=args.embryos, hsquared_broad=args.heritability_broad, 
            hsquared_snp=args.heritability_snp, rsquared=args.attenuated_rsquared, r_g=args.correlation_pleio, hsquared_pleio=args.heritability_pleio)
        
        # generate indices
        rand_index = get_random_index(args.embryos, args.n)
        max_index = get_max_pgs_index(values['ghat'])
        
        # get max and random EA values
        max_y = select_embryos_by_index(values['y'], max_index)
        rand_y = select_embryos_by_index(values['y'], rand_index)
        
        # get matching max and random pleio values
        max_pleio = select_embryos_by_index(values['pleio'], max_index)
        rand_pleio = select_embryos_by_index(values['pleio'], rand_index)
        
        # calculate difference for EA just to verify
        diffs_y = calc_phenotype_diffs(max_y, rand_y, args.sd_pheno)
        mean = diffs_y['selected_values'].mean()
        interval = diffs_y['selected_values'].std() * 1.96
        print('For parents of ' + str(args.ancestry) + ' ancestry, the within-family 95-percent prediction interval for the difference in ' + args.main_pheno +  \
            ' between random and selected embryos is ' + '%.2f' % mean + ' +/- ' +'%.2f' % interval + '.')
        
        # calculate pleiotropic liability
        df_liability = pd.DataFrame()
        df_liability['max'] = np.random.normal(loc=max_pleio['selected_values'], scale=(1-args.heritability_pleio)**0.5, size=int(args.n))
        df_liability['rand'] = np.random.normal(loc=rand_pleio['selected_values'], scale=(1-args.heritability_pleio)**0.5, size=int(args.n))
        
        # calculate whether they have phenotype
        df_binary = pd.DataFrame()
        df_binary['max'] = (df_liability['max'] > norm.ppf(1-args.prevalence_pleio)).astype(int)
        df_binary['rand'] = (df_liability['rand'] > norm.ppf(1-args.prevalence_pleio)).astype(int)
        
        # calculate prevalence
        pleio_prev_max = df_binary['max'].mean()*100
        pleio_prev_rand = df_binary['rand'].mean()*100
        
        # print results
        print('For parents of ' + str(args.ancestry) + ' ancestry, ' + args.pleiotropic_trait + ' has a within-family prevalence of ' + '%.2f' % pleio_prev_max + \
            '% in embryos selected using the polygenic score for ' + args.main_pheno + ' and a within-family prevalence of ' + '%.2f' % pleio_prev_rand + '% for randomly selected embryos.')
        print('End of analysis.')


