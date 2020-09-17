"""
Purpose: To simulate expected educational attainment gains from embryo selection between families.
Date: 10/09/2019
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
import argparse


def calc_between_family_values(n, no_embryos, hsquared_bfsnp, eur_bf_rsquared):
    """
    Purpose: to get the ghat_i and y_i for each family pair, where i={1,...,no_embryos}.

    Arguments:
        n: number of comparisons to make
        no_embryos: number of embryos for each comparison
        hsquared_bfsnp: naive "SNP heritability" (for between family comparisons)
        eur_bf_rsquared: R^2 (ancestry-specific)

    Returns:
        {'ghat':df_ghat, 'educ':df_y}. 
        Each dataframe has size (n x no_embryos) and holds polygenic scores and phenotype values (unscaled), respectively.
    """
    df_g = pd.DataFrame()
    df_ghat = pd.DataFrame()
    df_y = pd.DataFrame()

    # generate 1e6 values of ghat_i, y_i for each i. Note: g is used to calculate ghat and y.
    for i in range(no_embryos):
        df_g['g_{}'.format(i)] = np.random.normal(loc=0, scale=hsquared_bfsnp**0.5, size=int(n))
        df_ghat['ghat_{}'.format(i)] = np.random.normal(loc=df_g['g_{}'.format(i)], \
            scale=(hsquared_bfsnp**2/eur_bf_rsquared - hsquared_bfsnp)**0.5, size=int(n))
        df_y['y_{}'.format(i)] = np.random.normal(loc=df_g['g_{}'.format(i)], \
            scale=(1-hsquared_bfsnp)**0.5, size=int(n))
    
    return {'ghat':df_ghat, 'y':df_y}


def get_random_index(no_embryos, no_observations):
    """
    Purpose: to generate a list of random integers

    Arguments:
        no_embryos: max value of random integers generated, i.e. random integers will be generated in the range [0, no_embryos)
        no_observations: number of random integers to return

    Returns:
        random integers from the “discrete uniform” distribution of size no_observations and in the range of [0, no_embyros)
    """
    return np.random.randint(no_embryos, size=int(no_observations))


def get_max_pgs_index(df_ghat):
    """
    Purpose: to identify the column that holds the max polygenic score for each row

    Arguments:
        df_ghat: dataframe of polygenic scores, where each row is a new parent pair. 

    Returns:
        series of indices of the highest polygenic score per parent pair
    """
    return df_ghat.idxmax(axis=1).map(lambda x: int(x.lstrip('ghat_')))


def select_embryos_by_index(df, index_of_embryos):
    """
    Purpose: to select values of a dataframe using a valid set of indices

    Arguments: 
        df: dataframe from which you want to select values
        index_of_embryos: indices you are using to select values from df

    Returns: 
        dataframe of selected values
    """
    df_values = pd.DataFrame()
    relevant_values = []
    count = 0
    # get relevant phenotype values (and ghat) based on index
    for row in df.itertuples(index=False):
        relevant_values.append(row[index_of_embryos[count]])
        count += 1
    df_values['selected_values'] = relevant_values
    return df_values


def calc_phenotype_diffs(df_selected_embryos, df_random_embryos, sd_pheno):
    """
    Purpose: to calculate the difference in educational attainment of a randomly selected
    embryo vs. an embryo selected by highest PRS for educational attainment.

    Arguments:
        df_selected_embryos: dataframe of embryos selected by highest PRS
        df_random_embryos: dataframe of randomly selected embryos
        sd_pheno: standard deviation of education in a population

    Returns:
        dataframe of difference in education level (measured in years) between randomly selected embryos
        and those selected on highest PRS
    """
    return (df_selected_embryos - df_random_embryos)*sd_pheno


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
    parser.add_argument("--n", default=1e6, type=float, \
                        help="Number of parent pairs to simulate. Defaults to 1e6.")
    parser.add_argument("--embryos", default=10, type=int, \
                        help="Number of embryos from which to choose. Defaults to 10.")
    parser.add_argument("--hsquared_bf", default=0.2, type=float, \
                        help="Naive heritability, or between-family heritability. Defaults to 0.2.")
    parser.add_argument("--eur_bf_rsquared", default=0.1, type=float, \
                        help="Naive R^2, or between-family Rsquared for EUR ancestry. Defaults to 0.1.")
    parser.add_argument("--sd_pheno", default=3.2, type=float, \
                        help="Standard deviation of phenotype of interest. Defaults to 3.2 for years of education.")
    parser.add_argument("--scale_AMR", default=1.6, type=float, \
                        help="Factor to convert EUR R2 values to AMR. Defaults to 1.6.")
    parser.add_argument("--scale_EAS", default=2.0, type=float, \
                        help="Factor to convert EUR R2 values to EAS. Defaults to 2.0.")
    parser.add_argument("--scale_AFR", default=4.5, type=float, \
                        help="Factor to convert EUR R2 values to AFR. Defaults to 4.5.")
    return parser.parse_args()






def main():
    # import arguments
    args = process_arguments()

    # dictionary of rsquared (between family) values
    DICT_RSQUARED_BF = {'EUR':args.eur_bf_rsquared, 
        'AMR':args.eur_bf_rsquared/args.scale_AMR, 
        'EAS':args.eur_bf_rsquared/args.scale_EAS, 
        'AFR':args.eur_bf_rsquared/args.scale_AFR}

    for ancestry in DICT_RSQUARED_BF:
        # calculate values using ancestry-specific rsquared
        values = calc_between_family_values(args.n, args.embryos, args.hsquared_bf, DICT_RSQUARED_BF.get(ancestry))
        
        # generate indices
        rand_index = get_random_index(args.embryos, args.n)
        max_index = get_max_pgs_index(values['ghat'])
        
        # get max pheno values
        max_y = select_embryos_by_index(values['y'], max_index)
        
        # get random pheno values
        rand_y = select_embryos_by_index(values['y'], rand_index)
        
        # calculate difference
        diffs_y = calc_phenotype_diffs(max_y, rand_y, args.sd_pheno)
        mean = diffs_y['selected_values'].mean()
        interval = diffs_y['selected_values'].std() * 1.96
        print('For ' + str(ancestry) + ' ancestry, the between-family 95-percent prediction interval ' + \
            'for the phenotype of interest is ' + '%.2f' % mean + ' +/- ' +'%.2f' % interval + '.')

    pass


if __name__ == "__main__":
    main()

