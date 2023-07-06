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
from Droplet import Droplet

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
        volume_um3 = Droplet.calculate_volume(size=size_droplet_um, size_type=size_droplet_type)
        volume_microliter = Droplet.convert_um3_to_ul(volume_um3)
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
        efficiencies_dntp_type='square',
        efficiencies_access=None,
        efficiencies_denature=None,
        efficiencies_anneal=None,
        efficiencies_k=None,
        use_tensor=False,
        t_denature_cycle=30,
        t_anneal_cycle=62,
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
            if efficiencies_dntp_type == 'linear':
                dict_efficiencies_dntp = {k:counter_dnTP_at_start_of_cycle[k]/v for k, v in counter_dnTP_initial.items()}
            if efficiencies_dntp_type == 'square':
                dict_efficiencies_dntp = {k:(counter_dnTP_at_start_of_cycle[k]/v)**2 for k, v in counter_dnTP_initial.items()}
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
        ret['efficiencies_denature_curr'] = efficiencies_denature
        ret['efficiencies_anneal_curr'] = efficiencies_anneal
        ret['efficiencies_dntp_denature_curr'] = efficiencies_denature
        ret['efficiencies_k_curr'] = efficiencies_k

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
            counter_ATCG_combined_curr = counter_ATCG_combined_curr + counter_amplicon_ATCG_curr
            counters_amplicon_ATCG_curr.append(counter_amplicon_ATCG_curr)
        ret['counter_ATCG_combined_curr'] = counter_ATCG_combined_curr

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

    @staticmethod
    def run_droplet_PCR(
            amplicons,
            nums_template_initial=None,
            size_droplet_um=20,
            num_pcr_cycles:int=30,
            dnTP_mM:float=0.2,
            dna_ng_per_ul=10,
            fragment_size_avg:int=550, #bp
            primer_F_uM = 0.5,
            primer_R_uM = 0.5,
            t_denature_intial=95,
            t_denature_cycle=30,
            t_anneal_cycle=62,
            efficiencies_dntp=None,
            efficiencies_dntp_type='square',
            efficiencies_access=None,
            efficiencies_denature=None,
            efficiencies_anneal=None,
            efficiencies_k=None,
            library_type='dsDNA',
            use_tensor=False,
        ):
        """
        """
        ret = {}

        # process amplicons
        assert len(amplicons) > 0, 'please provide amplicon sequences'
        num_unique_amplicons = len(amplicons)
        lengths_amplicon = np.array(list(map(len, amplicons)))
        counters_amplicon_ATCG, counter_ATCG_combined = PCR.get_base_counts_from_sequences(amplicons)

        if not np.any(nums_template_initial):
            nums_template_initial = np.array([1]*num_unique_amplicons)
        
        nums_template_at_start_of_cycle = nums_template_initial
        ret['nums_template_initial'] = nums_template_initial


        # process droplet parameters
        volume_um3 = Droplet.calculate_volume(size=size_droplet_um, size_type='diameter')
        droplet_volume_ul = Droplet.convert_um3_to_ul(volume_um3)

        # process dnTPs
        num_dnTP = PCR.approximate_dntp_count_from_dntp_mM(dnTP_mM=dnTP_mM, volume_microliter=droplet_volume_ul)
        counter_dnTP_initial = Counter({'A':num_dnTP, 'T':num_dnTP, 'C':num_dnTP, 'G':num_dnTP})
        counter_dnTP_at_start_of_cycle = counter_dnTP_initial


        ret['nums_template_at_start_of_cycle'] = nums_template_at_start_of_cycle
        ret['counter_dnTP_at_start_of_cycle'] = counter_dnTP_at_start_of_cycle
        ret['counter_dnTP_initial'] = counter_dnTP_initial
        ret['amplicons'] = amplicons

        # process dna
        dna_ng_in_droplet = dna_ng_per_ul*droplet_volume_ul
        num_dnMP = PCR.approximate_nt_count_from_dna_ng(
                    dna_ng=dna_ng_in_droplet,
                    fragment_size_avg=fragment_size_avg,
                    library_type=library_type
                )

        # process primers

        # run cycles
        ncycle2ret = {}
        for n in range(num_pcr_cycles):
            ret_cycle = PCR.pcr_cycle(
                            nums_template_at_start_of_cycle=nums_template_at_start_of_cycle,
                            counter_dnTP_at_start_of_cycle=counter_dnTP_at_start_of_cycle,
                            counter_dnTP_initial=counter_dnTP_initial,
                            amplicons=amplicons,
                            t_denature_cycle=t_denature_cycle,
                            t_anneal_cycle=t_anneal_cycle,
                            efficiencies_dntp=efficiencies_dntp,
                            efficiencies_dntp_type=efficiencies_dntp_type,
                            efficiencies_access=efficiencies_access,
                            efficiencies_denature=efficiencies_denature,
                            efficiencies_anneal=efficiencies_anneal,
                            efficiencies_k=efficiencies_k,
                        )

            ncycle2ret[n] = ret_cycle
            nums_template_at_start_of_cycle = ret_cycle['nums_template_at_end_of_cycle']
            counter_dnTP_at_start_of_cycle = ret_cycle['counter_dnTP_at_end_of_cycle']

        ret['ncycle2ret'] = ncycle2ret

        return ret

    @staticmethod
    def cartprod(*arrays):
        """
        cartesian product
        """
        N = len(arrays)
        return np.transpose(np.meshgrid(*arrays, indexing='ij'),
                        np.roll(np.arange(N + 1), -1)).reshape(-1, N)

    @staticmethod
    def compartmentalize_bulk_to_drops(
            nums_template_amplicon,
            n_droplets,
            V_total_ul,
            Ws_dropids=None,
        ):
        """
        """
        ret = {}
        aid2dropids = {}
        dropid_to_aid2count = {}

        for aid, num_template_amplicon in enumerate(nums_template_amplicon):
            if np.any(Ws_dropids):
                a=list(range(n_droplets))
                aid2dropids[aid] = np.random.choice(a, size=num_template_amplicon, p=Ws_dropids)
            else:
                aid2dropids[aid] = np.random.randint(low=0, high=n_droplets, size=num_template_amplicon, dtype=int)
            for dropid in aid2dropids[aid]:
                if dropid not in dropid_to_aid2count:
                    dropid_to_aid2count[dropid] = {}
                if aid not in dropid_to_aid2count[dropid]:
                    dropid_to_aid2count[dropid][aid] = 1
                else:
                    dropid_to_aid2count[dropid][aid] += 1

        ret['aid2dropids'] = aid2dropids
        ret['dropid_to_aid2count'] = dropid_to_aid2count

        return ret


    @staticmethod
    def multiplex_dd_pcr(
            amplicons,
            nums_template_amplicon,
            V_i_bulk_ul=200,
            n_droplets=10000,
            determine_n_droplets_from_V_i=True,
            efficiency_bulk_to_droplet=0.8,
            size_droplet_um=20,
            size_type='diameter',
            func_droplet_size=None, #np.random.normal,
            scale_droplet_size=0,
            C_i_dna_ng_per_ul=10,
            C_i_dnTP_mM=0.2,
            C_i_primer_F_uM=0.5,
            C_i_primer_R_uM=0.5,
            C_i_polymerase_U=1, # to be updated
            library_type='dsDNA',
            fragment_size_avg=550,
            num_pcr_cycles=40,
            t_denature_intial=95,
            t_denature_cycle=30,
            t_anneal_cycle=62,
            efficiencies_dntp=None,
            efficiencies_dntp_type='square',
            efficiencies_access=None,
            efficiencies_denature=None,
            efficiencies_anneal=None,
            efficiencies_k=None,
        ):

        ret = {}


        # droplet size(s)
        V_droplet_um3 = Droplet.calculate_volume(size=size_droplet_um, size_type=size_type)
        V_droplet_ul = Droplet.convert_um3_to_ul(V_droplet_um3)

        # V_total
        if determine_n_droplets_from_V_i:
            n_droplets = int(V_i_bulk_ul*efficiency_bulk_to_droplet/V_droplet_ul)
        V_total_ul = V_droplet_ul * n_droplets
        ret['n_droplets'] = n_droplets
        ret['Ws_dropid'] = None

        if func_droplet_size:
            sizes_droplet = np.random.normal(size_droplet_um, scale_droplet_size, n_droplets)
            f = lambda x: Droplet.calculate_volume(x)
            Vs_droplet_um3 = f(sizes_droplet)
            f = lambda x: Droplet.convert_um3_to_ul(x)
            Vs_droplet_ul = f(Vs_droplet_um3)
            V_total_ul = sum(Vs_droplet_ul)
            ret['Vs_droplet_ul'] = Vs_droplet_ul
            Ws_dropid = np.array(Vs_droplet_ul) / V_total_ul
            ret['Ws_dropid'] = Ws_dropid         
        ret['V_total_ul'] = V_total_ul

        # amplicons
        assert len(amplicons) > 0, 'please provide at least one amplicon sequences'
        num_unique_amplicons = len(amplicons)
        lengths_amplicon = np.array(list(map(len, amplicons)))
        amplicon = np.array(amplicons)

        nums_template_amplicon = np.array(nums_template_amplicon)
        nums_template_amplicon = nums_template_amplicon * V_total_ul / V_i_bulk_ul
        nums_template_amplicon = nums_template_amplicon.astype('int')
        ret['nums_template_amplicon'] = nums_template_amplicon

        # dna
        mass_dna_ng = C_i_dna_ng_per_ul * V_total_ul
        if func_droplet_size:
            masses_dna_ng = mass_dna_ng * Vs_droplet_ul / V_total_ul

        # compartmentalize
        ret_bulk2drop = PCR.compartmentalize_bulk_to_drops(
                            nums_template_amplicon=nums_template_amplicon,
                            n_droplets=n_droplets,
                            V_total_ul=V_total_ul,
                            Ws_dropids=ret['Ws_dropid'],
                        )
        num_droplets_nonempty = len(ret_bulk2drop['dropid_to_aid2count'])
        num_droplets_empty = n_droplets - num_droplets_nonempty
        ret['num_droplets_nonempty'] = num_droplets_nonempty
        ret['num_droplets_empty'] = num_droplets_empty

        if func_droplet_size:
            sizes_droplet_nonemtpy = 'Vs_total_ul'[list(ret_bulk2drop['dropid_to_aid2count'].keys())]
        else:
            sizes_droplet_nonempty = np.array([V_droplet_ul]*len(ret_bulk2drop['dropid_to_aid2count']))
        
        combos_in_non_empty_droplets =  [tuple(k for k, v in d.items() for i in range(v)) 
                                for d in ret_bulk2drop['dropid_to_aid2count'].values()]
        combo2count = Counter(combos_in_non_empty_droplets)
        ret['combo2count'] = combo2count

        # run PCR
        combo2ret_pcr = {}
        combo2df = {}
        ls = []
        for k, v in combo2count.items():
            amplicons_droplet = [amplicons[i] for i in k]
            amplicon2count = Counter(amplicons_droplet)
            amplicons_droplet = list(amplicon2count.keys())
            nums_template_initial_droplet = list(amplicon2count.values())
            n_droplets_combo = v

            ret_pcr_droplet = PCR.run_droplet_PCR(
                                    amplicons=amplicons_droplet,
                                    nums_template_initial=nums_template_initial_droplet,
                                    size_droplet_um=size_droplet_um,
                                    num_pcr_cycles=num_pcr_cycles,
                                    dnTP_mM=C_i_dnTP_mM,
                                    dna_ng_per_ul=C_i_dna_ng_per_ul,
                                    fragment_size_avg=fragment_size_avg, #bp
                                    primer_F_uM = C_i_primer_F_uM,
                                    primer_R_uM = C_i_primer_R_uM,
                                    t_denature_intial=t_denature_intial,
                                    t_denature_cycle=t_denature_cycle,
                                    t_anneal_cycle=t_anneal_cycle,
                                    efficiencies_dntp=efficiencies_dntp,
                                    efficiencies_dntp_type=efficiencies_dntp_type,
                                    efficiencies_access=None,
                                    efficiencies_denature=None,
                                    efficiencies_anneal=None,
                                    efficiencies_k=efficiencies_k,
                                    library_type='dsDNA',
                                    use_tensor=False,
                                )
            combo2ret_pcr[k] = ret_pcr_droplet

            df_a_droplet = pd.DataFrame.from_dict(ret_pcr_droplet['ncycle2ret'], orient='index')

            dfcombo = df_a_droplet[['amplicons']].copy()
            dfcombo['nums_template_at_start_of_cycle'] = df_a_droplet['nums_template_at_start_of_cycle'].apply(lambda x: np.array(x)*n_droplets_combo)
            dfcombo['{}_amplicon2counter'.format(k)] = dfcombo.apply(lambda x: Counter(dict(zip(x['amplicons'], x['nums_template_at_start_of_cycle']))), axis=1)
            ls.append(dfcombo['{}_amplicon2counter'.format(k)])

            combo2df[k] = dfcombo
        
        df_pooled = pd.concat(ls, axis=1)
        df_pooled['all'] = df_pooled.apply(lambda x: sum(x, Counter()), axis=1)
        df_pooled['V_total_ul'] = V_total_ul
        ret['df_pooled'] = df_pooled

        ret['combo2ret_pcr'] = combo2ret_pcr
        ret['combo2df'] = combo2df

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
        

    