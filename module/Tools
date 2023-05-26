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


    @method
    def calculate_volume(size:float, input_type:str='diameter'):
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

    @method
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

    @method
    def convert_um3_to_ml(x):
        return x*math.pow(10,-15+3)

    @method
    def convert_um3_to_ul(x):
        return x*math.pow(10,-15+6)

    @method
    def convert_um3_to_nl(x):
        return x*math.pow(10,-15+9)

    @method
    def convert_um3_to_pl(x):
        return x*math.pow(10,-15+12)

    @method
    def calculate_num_cell_at_saturation(size:float,
                                         size_type:str='diameter',
                                         saturation_density=math.pow(10, 10)):
        
        
        """
        # saturation density at 10^9
        # convert_um3_to_ml(calculate_volume_sphere(30))*math.pow(10,9)
        d_droplet: radius of the droplet in um
        sat_cell: saturation density of the cell, default 10^10 cells/ml
        """
        if size_type == 'diameter':
            size = size / 2

        V= Droplet.calculate_volume(size, size_type=size_type)
        V = Droplet.convert_um3_to_ml(V)
        num_cells = V * saturation_density
        num_cells = int(num_cells)

        return num_cells