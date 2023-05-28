import os
import math
import flowkit as fk
import numpy as np
import pandas as pd
import random
from collections import Counter
import matplotlib.pyplot as plt
import itertools

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
        n_rounds:int=10000, 
        colname_f1:str='mCherry-A', 
        colname_strain:str='sid',
        colname_strainP:str='P_sampling_sid',
        colname_indexP:str='P_sampling_index',
        max_num_cells_at_saturation:int=300):
        """
        @param
        
        """
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
        df_sampled= df.sample(n=n_rounds, weights=df['P_sampling_index'])
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
            num_cells_encapsulated=1,
            max_num_cells_at_saturation:int=300):
        """
        """
  

        assert num_cells_encapsulated >= 0, print('num_cells_encapsulated shall be an positive interger')

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
        ls_dfs = []
        if num_cells_encapsulated == 0:
            empty_record = [[], [], ()]
            cols = ['indicies_padded', 'values_padded', 'sid']
            dfs = pd.DataFrame([empty_record]*n_rounds, columns=cols)
            ls_dfs.append(dfs)
            strainscombo2count = Counter({():n_rounds})

        for rs in range(num_cells_encapsulated):
            df_sampled = df.sample(n=n_rounds, weights=df['P_sampling_index'])
            ls_Series.append(df_sampled[colname_strain].to_list())
            S_strains = pd.DataFrame(ls_Series).apply(','.join).apply(lambda x: x.split(','))
            S_strains = S_strains.apply(sorted).apply(tuple)
            strainscombo2count = Counter(S_strains)

        # sample cells in each droplet
        for combo, n1 in strainscombo2count.items():
            if combo == (): break
            strain2count = Counter(combo)
            dfstrain = df[df[colname_strain].isin(combo)].copy()

            for strain, n2 in strain2count.items():
                dfstrain.loc[dfstrain[colname_strain] == strain, 'P_sampling_index']*=n2
                Psum = dfstrain['P_sampling_index'].sum()
                dfstrain.loc[:, 'P_sampling_index'] /= Psum

            # use sampling with replace if n > than dataset
            replace = False
            if n1*max_num_cells_at_saturation > dfstrain.shape[0]:
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
        ret['df']['num_cells_encapsulated'] = num_cells_encapsulated
        
        return ret
    
    @staticmethod
    def simsort(
        df, 
        n_rounds:int=10000,
        colname_f1:str='mCherry-A',
        colname_strain:str='sid',
        colname_strainP:str='P_sampling_sid',
        colname_indexP:str='P_sampling_index',
        size_droplet=20,
        size_type='diameter',
        func_droplet_size=np.random.normal,
        scale_droplet_size:float=0,
        size_left_curve_only=False,
        bins:int=100,
        num_cells_encapsulated:int=1,
        func_cells_encapsulated_per_droplet=np.random.poisson,
        cell_encapsulation_rate:float=0.1,
        discard_empty_droplets:bool=False,
        rng=np.random.default_rng(),
        figsize=(5,5)):

        """
        @param 
        """
  
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

        # varying droplet_size
        if func_droplet_size:
            if size_left_curve_only:
                ls_size_droplet = func_droplet_size(size_droplet, 
                                                    size_droplet*scale_droplet_size, 
                                                    n_rounds*3)
                ls_size_droplet = [sz for sz in ls_size_droplet if sz <= size_droplet]
                ls_size_droplet = ls_size_droplet[:n_rounds]
            else:
                ls_size_droplet = func_droplet_size(size_droplet, 
                                                    size_droplet*scale_droplet_size, 
                                                    n_rounds)
            ls_num_cells_at_saturation_in_droplet = \
                list(map(dp.calculate_num_cell_at_saturation, ls_size_droplet))
        else:
            ls_size_droplet = [size_droplet]*n_rounds
            ls_num_cells_at_saturation_in_droplet  = \
                [dp.calculate_num_cell_at_saturation(size_droplet)]*n_rounds
        
        # get maximum droplet_size
        max_num_cells_at_saturation = max(ls_num_cells_at_saturation_in_droplet)
        ret['max_num_cells_at_saturation'] = max_num_cells_at_saturation

        ls_volume_droplet_t0_um3= list(map(dp.calculate_volume, ls_size_droplet))
        ls_volume_droplet_t0_pl = list(map(dp.convert_um3_to_pl, ls_volume_droplet_t0_um3))

        # varying num of encapsulated cells in each droplet
        if func_cells_encapsulated_per_droplet:
            a_num_cells_encapsulated = func_cells_encapsulated_per_droplet(cell_encapsulation_rate, n_rounds)
            
            # whether to drop empty cells
            if discard_empty_droplets:
                P1 = cell_encapsulation_rate*math.exp(-cell_encapsulation_rate)
                a_num_cells_encapsulated = func_cells_encapsulated_per_droplet(
                                                cell_encapsulation_rate, int(n_rounds/(P1/2))
                                            )
                a_num_cells_encapsulated = a_num_cells_encapsulated[a_num_cells_encapsulated > 0][:n_rounds]

        else:
            a_num_cells_encapsulated = np.array([num_cells_encapsulated]*n_rounds)

        # get sorted droplets
        ls_dfs = []

        ret['cells_encapsulated2count'] = Counter(a_num_cells_encapsulated)
        for num_cells_encapsulated, count in ret['cells_encapsulated2count'].items():
        
            tmp_ret = DropletSorter.get_droplets_encapsulated_n_cells(
                            df=df,
                            n_rounds=count,
                            colname_f1=colname_f1, 
                            colname_strain=colname_strain,
                            colname_strainP=colname_strainP,
                            colname_indexP=colname_indexP,
                            num_cells_encapsulated=num_cells_encapsulated,
                            max_num_cells_at_saturation=max_num_cells_at_saturation)
            #print(num_cells_encapsulated, count, tmp_ret['df'].shape)
            ls_dfs.append(tmp_ret['df'])


        ret['df'] = pd.concat(ls_dfs, axis=0, ignore_index=True)
        ret['df']['size_droplet'] = ls_size_droplet
        ret['df']['v_droplet_pl_t0'] = ls_volume_droplet_t0_pl

        #set empty droplets reflect no cell growth
        ret['df']['num_cells_at_saturation_in_droplet'] = ls_num_cells_at_saturation_in_droplet
        ret['df']['num_cells_at_saturation_in_droplet'] = ret['df'].apply(lambda x:
            0 if x['num_cells_encapsulated'] == 0 else x['num_cells_at_saturation_in_droplet'],
            axis=1
            )
        
        ret['df']['sum_mCherry'] = ret['df'].apply(
            lambda x: sum(x['values_padded'][:x['num_cells_at_saturation_in_droplet']]),
            axis=1)
        
        #plot histogram
        cols=['sum_mCherry', 'size_droplet', 'num_cells_encapsulated', 'num_cells_at_saturation_in_droplet']
        ret['fig'], ret['ax'] = plt.subplots()
        ret['df'][cols].hist(bins=bins,figsize=figsize, ax=ret['ax'])

        return ret

    @staticmethod
    def merge_droplets(
        df_droplets, 
        df_original,
        colname_f1='mCherry-A',
        colname_strain='sid'):
    
        ret = {}

        S_indicies = df_droplets.apply(
            lambda x: x['indicies_padded'][:x['num_cells_at_saturation_in_droplet']], 
            axis=1)
        indicies = list(itertools.chain(*S_indicies))
        df = pd.DataFrame(indicies, columns=['idx'])
        df[colname_f1] = df['idx'].apply(lambda x: df_original.loc[x, colname_f1])
        df[colname_strain] = df['idx'].apply(lambda x: df_original.loc[x, colname_strain])

        return df.copy()