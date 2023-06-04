import os
import math
import flowkit as fk
import numpy as np
import pandas as pd
import random
from collections import Counter
import matplotlib.pyplot as plt
import itertools


class DNA(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        """
        super(DNA, self).__init__()
    
    @staticmethod
    def generate_random_dna_sequence(
        size=100,
        p=(0.25,0.25,0.25,0.25) # ATCG
    )->str:
        
        """
        Sample usage
        ---
        size=100,
        p=(0.25,0.25,0.25,0.25) # ATCG
        generate_random_dna_sequence(
            size=size,
            p=p
        )
        """

        seq = np.random.choice(
                a=('A','T','C','G'), 
                size=size,
                replace=True,
                p=p)
        
        return ''.join(seq)