import os
import math
import flowkit as fk
import numpy as np
import pandas as pd
import random
from collections import Counter
import matplotlib.pyplot as plt

#local libraries
from Data import FlowData as fld
from Tools import Droplet as dp

class DropletSorter(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        """
        super(DropletSorter, self).__init__()
    
    @staticmethod
    def get_droplets_encapsulated_a_cell(
      df, 
      n_rounds= 10000, 
      colname_f1:str='mCherry-A', 
      colname_strain:str='sid',
      colname_strainP:str='P_sampling_sid',
      colname_indexP:str='P_sampling_index',
      max_num_cells_at_saturation:int=300):
  
        ret = {}
                
        # get records for strain probabilities
        dftmp = df.groupby(colname_strain)[colname_strain].count().rename('count').reset_index()
        ret['strain2count'] = dict(zip(dftmp[colname_strain], dftmp['count']))
        if colname_strain in df.columns and colname_strainP in df.columns:
            dftmp = df[[colname_strain, colname_strainP]].drop_duplicates()
            ret['strain2P'] = dict(zip(dftmp[colname_strain], dftmp[colname_strainP]))
        else:
            tmp_ret = fld.assign_sampling_probability(df, colname_strain=colname_strain)
            df = tmp_ret['df']
            ret['strain2P'] = tmp_ret['strain2P']
            ret['strain2count'] = tmp_ret['strain2count']
        
        if colname_indexP not in df.columns:
            df['P_sampling_index'] = df.apply(
                    lambda x: ret['strain2P'][x[colname_strain]]/ret['strain2count'][x[colname_strain]],
                    axis=1)

        # Sample strain for each round
        df_sampled= dfc_wz.sample(n=n_rounds, weights=dfc_wz['P_sampling_index'])
        df_sampled_count = df_sampled.groupby('sid')['sid'].count().rename('count').reset_index()
        ret['strain2n'] = dict(zip(df_sampled_count['sid'], df_sampled_count['count']))

        ls_dfs = []
        for strain, n in ret['strain2n'].items():
            dfstrain = df[df['sid'] == strain]

            # sample # cells/strain for each droplet
            ls_index = list(dfstrain.index)

            # use sampling with replace if n > than dataset
            replace = False
            if n*max_num_cells_at_saturation > dfstrain.shape[0]:
                replace = True
            
            indices_in_droplets = np.random.choice(
                                    a=ls_index,
                                    size=(n, max_num_cells_at_saturation),
                                    replace=replace,
                                )
            
            dfs = pd.Series(list(indices_in_droplets)).to_frame('indicies_padded')
            dfs['values_padded'] = dfs['indicies_padded'].apply(
                lambda x: [df.loc[i, colname_f1] for i in x])
            dfs['sid'] = [(strain,)]*dfs.shape[0]
            ls_dfs.append(dfs)
            
        ret['df'] = pd.concat(ls_dfs).reset_index(drop=True)
        
        return ret
    
    @staticmethod
    def get_droplets_encapsulated_n_cells(
                df,
                n_rounds=10000,
                colname_f1:str='mCherry-A', 
                colname_strain:str='sid',
                colname_strainP:str='P_sampling_sid',
                colname_indexP:str='P_sampling_index',
                num_cells_encapsulated=0,
                max_num_cells_at_saturation:int=300):
        """
        """
  

        assert num_cells_encapsulated > 0, print('num_cells_encapsulated shall be an positive interger')

        ret = {}

        # get records for strain probabilities
        dftmp = df.groupby(colname_strain)[colname_strain].count().rename('count').reset_index()
        ret['strain2count'] = dict(zip(dftmp[colname_strain], dftmp['count']))
        if colname_strain in df.columns and colname_strainP in df.columns:
            dftmp = df[[colname_strain, colname_strainP]].drop_duplicates()
            ret['strain2P'] = dict(zip(dftmp[colname_strain], dftmp[colname_strainP]))
        else:
            tmp_ret = fld.assign_sampling_probability(df, colname_strain=colname_strain)
            df = tmp_ret['df']
            ret['strain2P'] = tmp_ret['strain2P']
            ret['strain2count'] = tmp_ret['strain2count']
        
        if colname_indexP not in df.columns:
            df['P_sampling_index'] = df.apply(
                    lambda x: ret['strain2P'][x[colname_strain]]/ret['strain2count'][x[colname_strain]],
                    axis=1)
            
        # sample strain combinations captured in each droplet
        ls_Series = []
        for rs in range(num_cells_encapsulated):
            df_sampled = dfc_wz.sample(n=n_rounds, weights=dfc_wz['P_sampling_index'])
            ls_Series.append(df_sampled[colname_strain].to_list())
        S_strains = pd.DataFrame(ls_Series).apply(','.join).apply(lambda x: x.split(','))
        S_strains = S_strains.apply(sorted).apply(tuple)
        strainscombo2count = Counter(S_strains)

        # sample cells in each droplet
        ls_dfs = []
        for combo, n1 in strainscombo2count.items():
            strain2count = Counter(k)
        dfstrain = df[df[colname_strain].isin(combo)].copy()

        for strain, n2 in strain2count.items():
            dfstrain.loc[dfstrain[colname_strain] == strain, 'P_sampling_index']*=n2
            Psum = dfstrain['P_sampling_index'].sum()
            dfstrain.loc[:, 'P_sampling_index'] /= Psum

            # use sampling with replace if n > than dataset
            replace = False
            if n_rounds*max_num_cells_at_saturation > dfstrain.shape[0]:
                replace = True

            ls_index = list(dfstrain.index)
            indices_in_droplets = np.random.choice(
                                        a=ls_index,
                                        size=(n1, max_num_cells_at_saturation),
                                        replace=replace,
                                    )
            dfs = pd.Series(list(indices_in_droplets)).to_frame('indicies_padded')
            dfs['values_padded'] = dfs['indicies_padded'].apply(
                    lambda x: [df.loc[i, colname_f1] for i in x])
            dfs['sid'] = [combo]*dfs.shape[0]
            ls_dfs.append(dfs.copy())

        ret['df'] = pd.concat(ls_dfs).reset_index(drop=True)
        
        return ret