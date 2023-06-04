import os
import math
import flowkit as fk
import numpy as np
import pandas as pd
import random
import collections
from collections import Counter
import matplotlib.pyplot as plt
import itertools


class PCR(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        """
        super(PCR, self).__init__()
    
    @staticmethod
    def get_dnTP_counts(
        size_droplet_um:float=20,
        size_droplet_type:str='diameter',
        dnTP_mM:float=0.2,
        dict_ATCG_mM:dict=None,
    )->collections.Counter:
        
        """
        Sample usage
        ----
        size_droplet_um=20
        size_droplet_type='diameter'
        dnTP_mM=0.2,
        dict_ATCG_mM=None,
        
        """
        volume_um3 = dp.calculate_volume(size=size_droplet_um, size_type=size_droplet_type)
        volume_microliter = dp.convert_um3_to_ul(volume_um3)
        num_dnTP_available = approximate_num_nt_from_dntp_mM(
        dnTP_mM=dnTP_mM, 
        volume_microliter=volume_microliter
        )
        Counter_ATCG_available = Counter({k:num_dnTP_available for k in ('A', 'T', 'C', 'G')})

        if isinstance(dict_ATCG_mM, dict):
            Counter_ATCG_available = {}
            for k, dnTP_mM in dict_ATCG_mM.items():
                Counter_ATCG_available[k] = approximate_num_nt_from_dntp_mM(
                                                dnTP_mM=dnTP_mM, 
                                                volume_microliter=volume_microliter
                                            )
            Counter_ATCG_available = Counter(Counter_ATCG_available)

        return Counter_ATCG_available

    
    @staticmethod
    def approximate_num_nt_from_dntp_mM(
        dnTP_mM,
        volume_microliter,
        molecular_weight_dATP:float=491.2,
        molecular_weight_dCTP:float=467.2,
        molecular_weight_dGTP:float=507.2,
        molecular_weight_dTTP:float=482.2,
        molecular_weight_dnTP:float=487.,
        )->float:

        """
        Default settings
        ----
        molecular_weight_dATP=491.2,
        molecular_weight_dCTP=467.2,
        molecular_weight_dGTP=507.2,
        molecular_weight_dTTP=482.2,
        molecular_weight_dnTP=487,

        Sample usage
        ----
        dnTP_mM=0.2
        volume_microliter=50
        approximate_num_nt_from_dntp_mM(
        dnTP_mM=dnTP_mM, 
        volume_microliter=volume_microliter
        )

        """
        dnTP_mmole_per_liter = dnTP_mM
        dnTP_mmole = dnTP_mmole_per_liter*volume_microliter*10e-6
        num_dnTP = dnTP_mmole * 10e-6 * scipy.constants.Avogadro

        return(num_dnTP)