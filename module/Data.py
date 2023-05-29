import pandas as pd
import numpy as np
from collections import Counter
import itertools
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

class FlowData(object):
  
    def __init__(self):
        super(FlowData, self).__init__()

    
    @staticmethod
    def combine_flowcytometry_records(
        session, 
        ls_sampleid, 
        ls_bin_columns=None, 
        source='raw'):
        
        ls_df = []
        for sid in ls_sampleid:
            df = session.get_sample(sid).as_dataframe(source=source)
            df[('sid', 'sid')] = sid
            ls_df.append(df)
        dfc = pd.concat(ls_df, axis=0)
        if ls_bin_columns:
            for colname in ls_bin_columns:
                dfc[colname] = 1
        dfc.columns = dfc.columns.droplevel(0)
        dfc = dfc.reset_index()
        return dfc.copy()


    @staticmethod
    def generate_abuduance(distribution:str='gamma', 
                        a:int=3, 
                        scale:int=5, 
                        size:int=10):
        
        rng=np.random.default_rng()
        if distribution == 'gamma':
            abundances = rng.gamma(a, scale=scale, size=size)
    
        abundances = abundances/abundances.sum()
        return abundances
    
    @staticmethod
    def assign_sampling_probability(df, 
                              colname_strain='sid', 
                              strain2P:dict=None,
                               ):
        ret = {}

        # assign probabilities to strain
        dftmp = df.groupby(colname_strain)[colname_strain].count().rename('count').reset_index()
        ret['strain2count'] = dict(zip(dftmp[colname_strain], dftmp['count']))

        cols = df.columns
        cols = [col for col in cols if col != 'P_sampling_{}'.format(colname_strain)]
        if not strain2P:
            dftmp['P_sampling_{}'.format(colname_strain)] = dftmp['count']/dftmp['count'].sum()
            df = df[cols].merge(dftmp[[colname_strain, 'P_sampling_{}'.format(colname_strain)]],
                        on=colname_strain, how='left')
            ret['strain2P'] = dict(zip(dftmp[colname_strain], dftmp['P_sampling_{}'.format(colname_strain)]))
        else:
            df['P_sampling_{}'.format(colname_strain)] = df[colname_strain].apply(lambda x: strain2P[x])
            ret['strain2P'] = strain2P


        # assign probabilities to each event
        df['P_sampling_index'] = df.apply(
            lambda x: ret['strain2P'][x[colname_strain]]/ret['strain2count'][x[colname_strain]],
            axis=1)
        ret['df'] = df.copy()

        return ret
    
    @staticmethod
    def get_strain_counts_in_gate(
        df,
        peak_x:float, 
        colname_f1:str='mCherry-A',
        colname_strain:str='sid',
        width=10):
  
        strains = df[(df[colname_f1] > peak_x - width/2) & (df[colname_f1] > peak_x + width/2)][colname_strain]

        if type(df[colname_strain].iloc[0]) is str:
            ret = Counter(list(itertools.chain(strains)))
        elif type(df[colname_strain].iloc[0]) is tuple:
            ret = Counter(list(itertools.chain(*strains)))

        return ret
    
    @staticmethod
    def identify_strains_at_peaks(
        df,
        colname_f1:str='mCherry-A',
        colname_strain:str='sid',
        pct_width:float=0.01,
        bins:int=100,
        distance:float=10,
        height:float=250,
        above_quantile:float=0.99):

        ret = {}
        width = (df[colname_f1].max() - df[colname_f1].min())*pct_width

        hist, bin_edges = np.histogram(df[colname_f1].to_numpy(), bins=bins)
        peaks, peak_heights = find_peaks(hist, distance=distance, height=height)
        bin_edges = bin_edges[1:]
        arr_peak_x = bin_edges[peaks]
    
        ret['peak_heights'] = peak_heights
        ret['peaks_x'] = arr_peak_x

        ls = []
        for peak_x in arr_peak_x:
            record = {}
            record['peak_x'] = peak_x
            record['width'] = width
            strain2count = FlowData.get_strain_counts_in_gate(df, peak_x, colname_f1=colname_f1, width=width)
            record['num_strains'] = len(strain2count)
            record['counter_strains'] = strain2count
            record['num_events'] = sum(strain2count.values())

            print('peak_x: {}, width: {}, # strains present: {}, # events: {}'.format(
                peak_x, width, len(strain2count), record['num_events']))
            ls.append(record)
        
        record = {}
        v_at_quantile = df[colname_f1].quantile(q=above_quantile)

        strains = df[df[colname_f1] > v_at_quantile][colname_strain]
        record['v_at_quantile'] = v_at_quantile
        record['width'] = df[colname_f1].max() - v_at_quantile
        if type(df[colname_strain].iloc[0]) is str:
            strain2count = Counter(list(itertools.chain(strains)))
        elif type(df[colname_strain].iloc[0]) is tuple:
            strain2count = Counter(list(itertools.chain(*strains)))

        record['counter_strains'] = strain2count
        record['num_strains'] = len(strain2count)
        record['num_events'] = sum(strain2count.values())
        print('quantile {}: {}, width: {}, # strains present: {}, # events: {}'.format(above_quantile, 
                v_at_quantile, record['width'], len(strain2count), record['num_events']))
        

        ret['fig'] = plt.plot(bin_edges, hist)
        plt.plot(bin_edges[peaks], hist[peaks], "x")
        plt.show()
        
        ret['df'] = pd.DataFrame(ls)

        return ret