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