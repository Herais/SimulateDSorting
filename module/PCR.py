import os
import math
import scipy
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
    def approximate_dntp_count_from_dntp_mM(
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
        approximate_dntp_count_from_dntp_mM(
            dnTP_mM=dnTP_mM, 
            volume_microliter=volume_microliter
        )

        """
        dnTP_mmole_per_liter = dnTP_mM
        dnTP_mmole = dnTP_mmole_per_liter*volume_microliter*10e-6
        num_dnTP = dnTP_mmole * 10e-6 * scipy.constants.Avogadro

        return(num_dnTP)
    
    @staticmethod
    def approximate_nt_count_from_dna_ng(
        dna_ng,
        fragment_size_avg:int=550, #bp
        library_type:str='dsDNA',
        molecular_weight_dAMP:float=331.2,
        molecular_weight_dCMP:float=307.2,
        molecular_weight_dGMP:float=347.2,
        molecular_weight_dTMP:float=322.2,
        molecular_weight_dnMP:float=339.5,
        molecular_weight_ssDNAends:float=79.,
        molecular_weight_dsDNAends:float=157.9
        )->float:

        """
        Default settings
        ----
        molecular_weight_dAMP=331.2,
        molecular_weight_dCMP=307.2,
        molecular_weight_dGMP=347.2,
        molecular_weight_dTMP=322.2,
        molecular_weight_dnMP=339.5,
        molecular_weight_ssDNAends=79,
        molecular_weight_dsDNAends=157.9

        Sample usage
        ----
        dna_ng=500
        fragment_size_avg=550 #bp
        library_type='dsDNA'
        approximate_nt_count_from_dna_ng(
            dna_ng=dna_ng,
            fragment_size_avg=fragment_size_avg,
            library_type=library_type,
            )

        """
        if library_type == 'dsDNA':
            dna_frag_nmole = dna_ng / (molecular_weight_dnMP*fragment_size_avg*2 + molecular_weight_dsDNAends)
            dnMP_nmole = dna_frag_nmole*fragment_size_avg*2
            
        elif library_type == 'ssDNA':
            dna_frag_nmole = dna_ng / (molecular_weight_dnMP*fragment_size_avg + molecular_weight_ssDNAends)
            dnMP_nmole = dna_frag_nmole*fragment_size_avg

        num_dnMP = dnMP_nmole * 10e-9 * scipy.constants.Avogadro
    
        return(num_dnMP)
    
    @staticmethod
    def calculate_array_lengths_sequences(
        array_amplicon_sequences
    )->np.array(str):
        """
        Sample usage
        ----
        s1 = 'GTCTTCTAAACCTATGGACGGTA'
        s2 = 'CCGAAGGGTTTTAAACTAGTGCTACGCTAGTC'
        s3 = 'GCTAGTGCATGCCGCGCAGCGTATCAGTATCTAC'
        array_amplicon_sequences = np.array([s1, s2, s3])

        calculate_array_lengths_sequences(
            array_amplicon_sequences=array_amplicon_sequences,
        )
        """
        array_lengths_amplicon = np.array(list(map(len, array_amplicon_sequences)))
        
        return array_lengths_amplicon

    @staticmethod
    def get_base_counts_from_sequences(
        sequences:list,
    ):
        """
        Sample usage:
        s1 = 'GTCTTCTAAACCTATGGACGGTA'
        s2 = 'CCGAAGGGTTTTAAACTAGTGCTACGCTAGTC'
        s3 = 'GCTAGTGCATGCCGCGCAGCGTATCAGTATCTAC'
        ls_counter, counter_combined = get_base_counts_from_sequences([s1,s2,s3])

        """
        ls_counters = list(map(Counter, sequences))
        counter_combined = sum(ls_counters, Counter())
        
        return ls_counters, counter_combined
    
    @staticmethod
    def calculate_efficiencies_dntp(
        Counter_ATCG_initial,
        Counter_ATCG_available,
    )->(Counter, float):

        """
        Counter_ATCG_initial = Counter({'A': 200, 'T': 200, 'G': 200, 'C': 200})
        Counter_ATCG_available = Counter({'A': 20, 'T': 30, 'G': 15, 'C': 50})

        calculate_efficiencies_dntp(
            Counter_ATCG_initial=Counter_ATCG_initial,
            Counter_ATCG_available=Counter_ATCG_available,
        )
        """

        assert all(np.array(list(Counter_ATCG_available.values())) > 0), print(
            'There is insufficient dntp to complete the current cycle of PCR,'
        )

        dict_ATGC_efficiencies = {k:Counter_ATCG_available[k]/v for k, v in Counter_ATCG_initial.items()}
        efficiency_dntp = min(dict_ATGC_efficiencies.values())

        return dict_ATGC_efficiencies, efficiency_dntp
    
    @staticmethod
    def convert_primer_uM_to_count(
        primer_uM:float=0.5,
        volume_ul:float=4.188790204786391e-06,
    )->float:

        """
        Sample usage
        ----
        primer_F_uM = 0.5
        primer_R_uM = 0.5
        volume_um3 = Droplet.calculate_volume(size=20)
        droplet_volume_ul=dp.convert_um3_to_ul(volume_um3)

        count_primer_F = convert_primer_uM_to_count(primer_uM=primer_F_uM, volume_ul=droplet_volume_ul)
        count_primer_R = convert_primer_uM_to_count(primer_uM=primer_R_uM, volume_ul=droplet_volume_ul)
        """
        primer_mole_per_liter = primer_uM * 1e-6
        primer_mole =  primer_mole_per_liter * volume_ul * 1e-6
        count_primers = primer_mole * scipy.constants.Avogadro

        return count_primers



class Polymerase(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        """
        super(Polymerase, self).__init__()


    def Taq(self):
        self.synthesize_speed_min_per_kb = 1
    
    def Phusion(self):
        self.synthesize_speed_min_per_kb = 2
        

    