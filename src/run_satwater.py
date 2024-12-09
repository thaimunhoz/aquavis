from src.satwater import build_satwater

def run_satwater(select_sat, tile, period_ini, period_end, output_dir):

    '''
    Run the SatWater processing chain.
    Input:
        select_sat (str): Satellite to process (landsat or sentinel)
        tile (str): Tile name
        period_ini (str): Initial date in the format 'YYYYMMDD'
        period_end (str): End date in the format 'YYYYMMDD'
        output_dir (str): Output directory
    Output:
        None
    '''

    tiles_sentinel = r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\MGRS_tiles.shp'
    tiles_landsat = r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\ms_landsat.shp'

    n_cores=12

    params = {}

    params['aux_info'] = {}
    params['aux_info']['period'] = (period_ini, period_end)
    params['aux_info']['n_cores'] = n_cores
    params['aux_info']['sat_name'] = select_sat

    params['sentinel'] = {}
    params['sentinel']['input_dir'] = r'Z:\dbcenter\images\sentinel\scenes\level_toa'
    params['sentinel']['tiles_shp'] = tiles_sentinel

    params['landsat'] = {}
    params['landsat']['input_dir'] = r'Z:\dbcenter\images\landsat\scenes\level_toa'
    params['landsat']['generation'] = 'L89'
    params['landsat']['tiles_shp'] = tiles_landsat

    params[select_sat]['tiles'] = tile
    params['output_dir'] = output_dir

    SatWater_i = build_satwater.SatWater(select_sat, params)

    # 1. Run atmospheric correction
    SatWater_i.run_atmcor()

    # 2. Run tiling - Landsat images are reprojected, resampled, and clipped to the sentinel tile (MGRS)
    if select_sat == 'landsat':
        SatWater_i.run_tiling()

    # 3. Run resampling - Sentinel images are resampled to the Landsat spatial resolution (30m)
    if select_sat == 'sentinel':
        SatWater_i.run_resample()

    # 4. Create the final HLS synthetic image
    SatWater_i.run_hlswater()

    # 5. Save true color composition
    SatWater_i.run_plot()

    print('Process finished.')