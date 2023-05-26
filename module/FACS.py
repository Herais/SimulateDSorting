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

class Flowcytometry(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(Flowcytometry, self).__init__()

    @staticmethod
    def analyze(df,
                n_rounds:int=10000,
                colname_f1:str='mCherry-A', 
                colname_strain:str='sid',
                colname_strainP:str='P_sampling_sid',
                colname_indexP:str='P_sampling_index',
                replace=False,
                bins=200,
                figsize=(5,5),
                ):
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


        ls_v = []

        ret['fig'], ret['ax'] = plt.subplots()

        idx =  df[colname_f1].sample(n_rounds, replace=replace, weights=df['P_sampling_index']).index
        ret['df'] = df.loc[idx,:].copy()
        ret['df'][colname_f1].hist(bins=bins,figsize=figsize, ax=ret['ax'])
        return ret