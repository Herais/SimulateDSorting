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

# import local libraries
from Dsort import DropletSorter

class PCR(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        """
        super(PCR, self).__init__()

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

        primer_F_copies = convert_primer_uM_to_count(primer_uM=primer_F_uM, volume_ul=droplet_volume_ul)
        primer_R_copies = convert_primer_uM_to_count(primer_uM=primer_R_uM, volume_ul=droplet_volume_ul)
        """
        primer_mole_per_liter = primer_uM * 1e-6
        primer_mole =  primer_mole_per_liter * volume_ul * 1e-6
        count_primers = primer_mole * scipy.constants.Avogadro

        return count_primers

    @staticmethod
    def calculate_efficiencies_dntp_for_amplicons(
            ls_counter_amplicon_ATCG,
            dict_efficiencies_dnTP,
        ):
        """
        Sample usage
        ----
        calculate_efficiencies_dntp_for_amplicons(ls_counter_amplicon_ATCG, dict_efficiencies_dntp)
        """

        ls = []
        for counter_amplicon_ATCG in ls_counter_amplicon_ATCG:
            total_dnTP_count = sum(counter_amplicon_ATCG.values())
            efficiency_amplicon_bases = np.array([counter_amplicon_ATCG[k]*v for k, v in dict_efficiencies_dnTP.items()])
            efficiency_amplicon = sum(efficiency_amplicon_bases)/total_dnTP_count
            ls.append(efficiency_amplicon)

        return np.array(ls)

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
        volume_um3 = DropletSorter.calculate_volume(size=size_droplet_um, size_type=size_droplet_type)
        volume_microliter = DropletSorter.convert_um3_to_ul(volume_um3)
        num_dnTP_available = PCR.approximate_num_nt_from_dntp_mM(
        dnTP_mM=dnTP_mM, 
        volume_microliter=volume_microliter
        )
        Counter_ATCG_available = Counter({k:num_dnTP_available for k in ('A', 'T', 'C', 'G')})

        if isinstance(dict_ATCG_mM, dict):
            Counter_ATCG_available = {}
            for k, dnTP_mM in dict_ATCG_mM.items():
                Counter_ATCG_available[k] = PCR.approximate_num_nt_from_dntp_mM(
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
        num_dnTP = dnTP_mmole * 1e-3 * scipy.constants.Avogadro

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
        molecular_weight_dAMP=331.2, # Dalton = 1g/mole = 1ng/nmole
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

        count_dnMP = dnMP_nmole * 1e-9 * scipy.constants.Avogadro
    
        return(count_dnMP)
    
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
    ):

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
    
    @staticmethod
    def pcr_cycle(
        nums_template_at_start_of_cycle:list,
        counter_dnTP_at_start_of_cycle:collections.Counter,
        counter_dnTP_initial:collections.Counter,
        amplicons,
        efficiencies_dntp=None,
        efficiencies_access=None,
        efficiencies_denature=None,
        efficiencies_anneal=None,
        efficiencies_k=None,
        use_tensor=False,
        T_denature=95,
        t_denature=30,
        T_anneal=62,
        t_anneal=15,
    )->dict:

        """
        """

        ret = {}
        ret['nums_template_at_start_of_cycle'] = nums_template_at_start_of_cycle
        ret['counter_dnTP_at_start_of_cycle'] = counter_dnTP_at_start_of_cycle
        ret['counter_dnTP_initial'] = counter_dnTP_initial
        ret['amplicons'] = amplicons


        nums_template_at_start_of_cycle = np.array(nums_template_at_start_of_cycle)

        assert len(amplicons) > 0, 'please provide amplicon sequences'
        num_unique_amplicons = len(amplicons)
        lengths_amplicon = np.array(list(map(len, amplicons)))

        counters_amplicon_ATCG, counter_ATCG_combined = PCR.get_base_counts_from_sequences(amplicons)


        if not np.any(efficiencies_dntp):
            assert len(counter_dnTP_at_start_of_cycle) == 4, 'At least one of the dnTPs has depleted'
            dict_efficiencies_dntp = {k:counter_dnTP_at_start_of_cycle[k]/v for k, v in counter_dnTP_initial.items()}
            efficiencies_dntp = PCR.calculate_efficiencies_dntp_for_amplicons(counters_amplicon_ATCG, dict_efficiencies_dntp)

        if not np.any(efficiencies_access):
            efficiencies_access = np.array([1.]*num_unique_amplicons)

        if not np.any(efficiencies_denature):
            efficiencies_denature = np.array([1.]*num_unique_amplicons)

        if not np.any(efficiencies_anneal):
            efficiencies_anneal = np.array([1.]*num_unique_amplicons)

        if not np.any(efficiencies_k):
            efficiencies_k = np.array([1.]*num_unique_amplicons)

        ret['efficiencies_dntp_curr'] = efficiencies_dntp
        ret['efficiencies_access_curr'] = efficiencies_access
        ret['efficiencies_dntp_denature_curr'] = efficiencies_denature
        ret['efficiencies_anneal_curr'] = efficiencies_anneal
        ret['efficiencies_dntp_denature_curr'] = efficiencies_denature

        # load to pytorch tensor
        if use_tensor:
            efficiencies_dntp = torch.tensor(efficiencies_dntp)
            efficiencies_denature = torch.tensor(efficiencies_denature)
            efficiencies_anneal= torch.tensor(efficiencies_anneal)
            efficiencies_k = torch.tensor(efficiencies_k)


        # calculate PCR output at the end of the current cycle
        nums_template_produced_curr = nums_template_at_start_of_cycle* \
                                    efficiencies_dntp* \
                                    efficiencies_access* \
                                    efficiencies_denature* \
                                    efficiencies_anneal* \
                                    efficiencies_k

        # ATCG consumed in current cycle
        counters_amplicon_ATCG_curr = []
        counter_ATCG_combined_curr = Counter()
        for counter_amplicon_ATCG, num_template in zip(counters_amplicon_ATCG, nums_template_produced_curr):
            counter_amplicon_ATCG_curr = Counter({k:v*num_template for k, v in counter_amplicon_ATCG.items()})
            counter_ATCG_combined = counter_ATCG_combined + counter_amplicon_ATCG_curr
            counters_amplicon_ATCG_curr.append(counter_amplicon_ATCG_curr)

        counter_dnTP_at_end_of_cycle = counter_dnTP_at_start_of_cycle - counter_ATCG_combined_curr
        num_dnTP_left = len(counter_dnTP_at_end_of_cycle)
        if num_dnTP_left < 4:
            ret['nums_template_at_end_of_cycle'] = nums_template_at_start_of_cycle
            ret['counter_dnTP_at_end_of_cycle'] = counter_dnTP_at_start_of_cycle
        else:
            nums_template_at_end_of_cycle = nums_template_at_start_of_cycle + nums_template_produced_curr
            ret['nums_template_at_end_of_cycle'] = nums_template_at_end_of_cycle
            ret['counter_dnTP_at_end_of_cycle'] = counter_dnTP_at_end_of_cycle
        
        
        return ret



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
        

    