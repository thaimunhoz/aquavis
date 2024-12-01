import os
from satwater import tiling
from satwater import hlswater
from satwater import visualization
import satwater.atmcor.atm6s as atmcor

class SatWater(object):

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

        tiling.tiles.run(self.params)

    def run_resample(self):

        tiling.resample.run(self.params)

    def run_hlswater(self):
        hlswater.generate_hls.run(self.params)

    def run_plot(self):
        visualization.plot_images.run(self.params)