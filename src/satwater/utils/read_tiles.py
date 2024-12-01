

import geopandas

def readshp(path_shp, select_sat='landsat'):

    feature = geopandas.read_file(path_shp)

    if select_sat=='sentinel':
        return feature['Name'].to_list()
    else:
        pathrow = [f"{i[0]:03d}_{i[1]:03d}" for i in zip(feature['PATH'].to_list(), feature['ROW'].to_list())]
        return pathrow