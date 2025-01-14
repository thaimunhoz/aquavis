from src.run_satwater import run_satwater

# Define parameters
select_sat = 'sentinel'
tile = '21HVB'
period_ini = '20230101'
period_end = '20230105'
output_dir = '/ddnlus/scratch/r3693/hls_water/HLS_DATASET/sentinel'

import sys
sys.path.append('/ddnlus/r3693/hls_water/scripts/hls_water/src')

# Run the function
run_satwater(select_sat, tile, period_ini, period_end, output_dir)
