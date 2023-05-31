import math

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
    def calculate_n_droplets_from_bulk(
        V_i:float=500, #ul
        size_droplet:float=20, #um
        efficiency:float=0.6,
        size_type:str='diameter')->int:

        volume_avg_droplet_um3 = dp.calculate_volume(
        size=size_droplet,
        size_type=size_type)

        volume_avg_droplet_ul = dp.convert_um3_to_ul(volume_avg_droplet_um3)
        n_droplets = int(V_i/volume_avg_droplet_ul*efficiency)

        return n_droplets