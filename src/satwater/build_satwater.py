import os
from src.satwater.atmcor import atm6s as atmcor
from src.satwater.tiling import tiles as tiling
from src.satwater.tiling import resample as resample
from .hlswater import generate_hls as hlswater
from .visualization import plot_images as visualization

class SatWater(object):

    '''
    Class to run the entire process of atmospheric correction, tiling, resampling, HLS water generation, and plotting.
    '''

    def __init__(self, select_sat=None, params=None):

        self.select_sat = select_sat
        self.params = params

    def _createdir(self):
        if not os.path.exists(self.params['output_dir']):
            os.makedirs(self.params['output_dir'])

    def run_atmcor(self):

        self._createdir()
        atmcor.run(self.select_sat, self.params)

    def run_tiling(self):

        tiling.run(self.select_sat, self.params)

    def run_resample(self):

        resample.run(self.params)

    def run_hlswater(self):
        hlswater.run(self.params)

    def run_plot(self):
        visualization.run(self.params)