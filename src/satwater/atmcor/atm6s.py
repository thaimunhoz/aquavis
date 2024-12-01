import os
import datetime
import satwater.atmcor.gceratmos_sentinel.run_gceratmos as gceratmos_sentinel

from multiprocessing import Pool

def checkdaterange(all_imgs, dates_period, select_sat='landsat'):

    """Check if images are within the given date range
     Args:
        all_imgs (list): List of images
        dates_period (list): List of start and end dates in the format 'YYYYMMDD'
        select_sat (str): Satellite type (default: 'landsat')
     Returns:
        list: List of images within the given date range
    """
    all_imgs_within = []
    for img in all_imgs:

        # Get the date string from the image name
        if select_sat == 'sentinel':
            date_str = os.path.basename(img).split('_')[2].split('T')[0]
        else:
            date_str = os.path.basename(img).split('_')[3]

        # Convert the date string to datetime object
        date_target = datetime.datetime.strptime(date_str, '%Y%m%d')
        start_date = datetime.datetime.strptime(dates_period[0], '%Y%m%d')
        end_date = datetime.datetime.strptime(dates_period[1], '%Y%m%d')

        # Check if the date is within the given range
        if start_date <= date_target <= end_date:
            all_imgs_within.append(img)

    return all_imgs_within

def run(select_sat, params):

    """Runs a given set of parameters to initiate a selection process for satellite data.
    The 'select_sat' parameter defines the satellite to select (defaults to 'landsat')"""

    tiles = params[select_sat]['tiles']

    for tile in tiles:

        # Get the image path filenames within period
        if select_sat == 'landsat':

            # Landsat
            src_dir_tile = fr'{params[select_sat]["input_dir"]}\{tile}\{params[select_sat]["generation"]}'
            all_imgs = [fr'{src_dir_tile}\{i}' for i in os.listdir(src_dir_tile)]

        else:

            #sentinel
            src_dir_tile = fr'{params[select_sat]["input_dir"]}\{tile}'
            all_imgs = [fr'{src_dir_tile}\{i}\{os.listdir(os.path.join(src_dir_tile, i))[0]}' for i in os.listdir(src_dir_tile)]

        # Check the period and target date
        all_imgs = checkdaterange(all_imgs, params['aux_info']['period'], select_sat=select_sat)  # check the image date within the period

        if not all_imgs:

            print(f'No images for {tile} in the period')

            return

        # output location
        output_rrs = [fr'{params["output_dir"]}\atmcor\{select_sat}\{tile}\{os.path.basename(i).split(".")[0]}' for i in all_imgs]
        select_sat_list = [select_sat] * len(all_imgs)

        args = zip(all_imgs, output_rrs, select_sat_list)

        #for i in range(len(all_imgs)):
        #    gceratmos_sentinel.run_gceratmos(all_imgs[i], output_rrs[i], select_sat)

        with Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap(gceratmos_sentinel.run_gceratmos, args)
            print(results)