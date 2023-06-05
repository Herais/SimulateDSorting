import math
import torch
import numpy as np
import pandas as pd

class Droplet(object):
 
    def __init__(self, size:int=20):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(Droplet, self).__init__()
        self.size = size
        self.volume = self.calculate_volume(size)
        self.surface_area = self.calculate_surface_area(size)


    @staticmethod
    def calculate_volume(size:float, size_type:str='diameter'):
        """
        Parameters
        ----------
        size : float
            radiu
        
        input_type = radius | diameter

        Returns
        -------
        v : float
            volume of the droplet
        """

        if size_type == 'diameter':
            size = size / 2

        v = 4/3*math.pi*size*size*size
        return v

    @staticmethod
    def calculate_surface_area(size:float, size_type:str='diameter'):
        """
        Parameters
        ----------
        size : float
            radiu
        
        input_type = radius | diameter

        Returns
        -------
        SA : float
            surface area of the droplet
        """
        if size_type == 'diameter':
            size = size / 2
        
        a = 4*math.pi*size*size
    
        return a

    @staticmethod
    def convert_um3_to_ml(x):
        return x*math.pow(10,-15+3)

    @staticmethod
    def convert_um3_to_ul(x):
        return x*math.pow(10,-15+6)

    @staticmethod
    def convert_um3_to_nl(x):
        return x*math.pow(10,-15+9)

    @staticmethod
    def convert_um3_to_pl(x):
        return x*math.pow(10,-15+12)

    @staticmethod
    def calculate_num_cell_at_saturation(size:float,
                                         size_type:str='diameter',
                                         saturation_density=math.pow(10, 10)):
        
        
        """
        # saturation density at 10^9
        # convert_um3_to_ml(calculate_volume_sphere(30))*math.pow(10,9)
        d_droplet: radius of the droplet in um
        sat_cell: saturation density of the cell, default 10^10 cells/ml
        """

        V= Droplet.calculate_volume(size, size_type=size_type)
        V = Droplet.convert_um3_to_ml(V)
        num_cells = V * saturation_density
        num_cells = int(num_cells)

        return num_cells

    @staticmethod
    def generate_droplets_with_variable_sizes(
        n_droplets=100000,
        size_droplet:float=20., #um
        size_type:str='diameter',
        func_droplet_size=torch.normal,
        std_droplet_size=0.1)->torch.tensor:
    
        ret = {}
        if func_droplet_size:
            sizes_droplet = torch.normal(
                mean=size_droplet, 
                std=torch.arange(std_droplet_size, n_droplets),
                )
        else:
            sizes_droplet = torch.tensor([size_droplet]*n_droplets)
        
        return sizes_droplet

    @staticmethod
    def calculate_n_droplets_from_bulk(
        V_i:float=500, #ul
        size_droplet:float=20, #um
        efficiency_encapsulation:float=0.8,
        size_type:str='diameter'
    )->int:

        """
        Sample Usage
        ----
        V_i=500 #ul
        size_droplet=20 #um
        efficiency_encapsulation=0.8
        size_type='diameter'

        calculate_n_droplets_from_bulk(
            V_i=V_i,
            size_droplet=size_droplet,
            efficiency_encapsulation=efficiency_encapsulation,
            size_type=size_type,
        )
        """

        volume_avg_droplet_um3 = dp.calculate_volume(
                                    size=size_droplet,
                                    size_type=size_type)

        volume_avg_droplet_ul = dp.convert_um3_to_ul(volume_avg_droplet_um3)
        n_droplets = int(V_i/volume_avg_droplet_ul*efficiency_encapsulation)

        return n_droplets
    
    @staticmethod
    def generate_droplets_with_variable_size(
        n_droplets:int=100000,
        size_droplet:float=20., #um
        size_type:str='diameter',
        func_droplet_size=torch.normal,
        std_droplet_size:float=0.1,
        return_numpy:bool=False,
    )->torch.tensor:
        """
        Sample Usage
        ----
        n_droplets=100000
        size_droplet=20. #um
        size_type='diameter'
        func_droplet_size=torch.normal
        std_droplet_size=0.1

        generate_droplets_with_variable_size(
            n_droplets=n_droplets,
            size_droplet=size_droplet,
            size_type=size_type,
            func_droplet_size=func_droplet_size,
            std_droplet_size=std_droplet_size
        )

        """
        ret = {}
        if func_droplet_size:
            sizes_droplet = torch.normal(
                mean=size_droplet, 
                std=torch.arange(std_droplet_size, n_droplets),
                )
        else:
            sizes_droplet = torch.tensor([size_droplet]*n_droplets)
        
        return sizes_droplet