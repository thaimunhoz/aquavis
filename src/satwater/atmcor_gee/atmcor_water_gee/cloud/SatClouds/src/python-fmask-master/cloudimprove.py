
'''
I did install https://visualstudio.microsoft.com/visual-cpp-build-tools/

and I added the path to fmask folder

'''

import os
import sys
import fmask

# Pattern to identify the SAFE directory of Sentinel 2 images
safe_dir_base = 'S2*_MSIL1C_*.SAFE'
# Base to save the cloud bands
base_cloud = '_cloud_mask.tif'
# Information regarding the cloud mask algorithm
shadow_buffer_distance = 300
cloud_buffer_distance = 600
pixel_size = 60
# Path where the conda environment was created
conda_env = r"/Users/rejanesp/anaconda3/envs/fmask_env/bin/python.exe"
# Path where the FMASK algorithm is located
fmask_algorithm = r"/Volumes/rsp/FMASK/python-fmask-master/bin/fmask_sentinel2Stacked.py"
path_safedir =r"/Users/rejanesp/Downloads/S2B_MSIL1C_20220901T161829_N0400_R040_T17TLH_20220901T182839.SAFE 2"
dir_output_date = r"/Volumes/rsp/FMASK/OUT"
temp_safe = r"/Volumes/rsp/FMASK/OUT"
path_file_cloud = fr"/Volumes/rsp/FMASK/OUT/cloud.tif"

# Run the FMASK algorithm
cmd = rf"{conda_env} {fmask_algorithm} --shadowbufferdistance {shadow_buffer_distance} --cloudbufferdistance " \
      rf"{cloud_buffer_distance} --pixsize {pixel_size} -e {temp_safe} -o {path_file_cloud} --safedir {path_safedir} --tempdir {dir_output_date}"

os.system(cmd)
#










#
#
#
#
#
# # def run_fmask(safename, temp_safe):
# #     import os
# #     import sys
# #     sys.path.append(r'Z:\projects\code_center\advance\fmask\python-fmask-0.5.8\fmask')
# #     sys.path.append(r'Z:\projects\code_center\advance\fmask\python-fmask-0.5.8\fmask\cmdline')
# #     import fmask
# #     print(f"Starting Fmask with {safename}")
# #     path_safedir1 = os.path.join(temp_safe, safename)
# #     nm = safename.replace(".SAFE", "")
# #     path_cloud = os.path.join(temp_safe, f"cloud_{nm}.tif")
# #
# #     path_safedir1=r'Y:\datasets\images\sentinel\level_toa\15SXR\20151121\S2A_MSIL1C_20151121T165612_N0204_R026_T15SXR_20151121T165611.SAFE'
# #     path_cloud = fr"C:\Users\vsm71\Downloads\cloud.tif"
# #     temp_safe = r'C:\Users\vsm71\Downloads\temp'
# #     cmd = rf'C:\anaconda3\envs\testenv\python.exe Z:\projects\code_center\advance\fmask\python-fmask-pythonfmask-0.5.8\bin\fmask_sentinel2Stacked.py --shadowbufferdistance 600 --cloudbufferdistance 300 --pixsize 180 -e {temp_safe} -o {path_cloud} --safedir {path_safedir1}'
# #     os.system(cmd)
#
#
#
#
#
#
#
# import sys, shutil
# sys.path.append('/mnt/research/sat.6/vitorsm/Sentinel_burn/scripts/Preprocessing/utils')
# import sentinelcloud
# import intersectraster
#
# def export_tif(output_path, burn_p, ds):
#
#     """Create GeoTiff image with 1 band
#
#             Parameters
#             ----------
#             output_path : str
#                 filename
#             burn_p : 2-D numpy array
#                 Raster with 1 band
#             ds : gdal object
#                 This object has all properties (geotransform etc)
#
#             Returns
#             -------
#             None
#
#             """
#
#     import gdal
#     # create image driver
#     driver = gdal.GetDriverByName('GTiff')
#     # create destination for label file
#     file = driver.Create(output_path,
#                          burn_p.shape[1],
#                          burn_p.shape[0],
#                          1,
#                          gdal.GDT_Byte, ['COMPRESS=LZW'])
#     file.SetGeoTransform(ds.GetGeoTransform())
#     file.SetProjection(ds.GetProjection())
#     #file.GetRasterBand(1).SetNoDataValue(0)
#     # write label file
#     file.GetRasterBand(1).WriteArray(burn_p)
#     file.FlushCache()  ##saves to disk!!
#     file = None
#     ds = None
#     driver = None
#
# def check_exist(outfolder):
#     if os.path.exists(outfolder):
#         return True
#
# def download_images(base_cloude, MGRS_TILE, data_inicial, data_final, outfolder):
#
#     import os
#     #import pandas as pd
#     # Filtering and pickle
#     #base_cloude = pd.read_csv('/mnt/ufs18/rs-sat.6/vitorsm/Sentinel_burn/scripts/download/indexL1.csv', encoding='utf-8')
#     #data_inicial = '2019-01-01'
#     #data_final = '2021-03-30'
#     #selecao = base_cloude[(base_cloude['SENSING_TIME'] >= data_inicial) & (base_cloude['SENSING_TIME'] < data_final)]
#     #selecao.to_pickle('/mnt/ufs18/rs-sat.6/vitorsm/Sentinel_burn/scripts/download/indexL1.pkl')  # where to save it, usually as a .pkl
#
#     #base_cloude = pd.read_pickle('/mnt/ufs18/rs-sat.6/vitorsm/Sentinel_burn/scripts/download/indexL1.pkl')
#
#     cobertura_nuvem = 100
#
#     print('Consultando pela data ...' + data_inicial + "...")
#     selecao = base_cloude[(base_cloude['SENSING_TIME'] >= data_inicial) & (base_cloude['SENSING_TIME'] < data_final)]
#     scenes = selecao[(selecao.MGRS_TILE == MGRS_TILE) & (selecao.CLOUD_COVER < cobertura_nuvem)]
#     scenes = scenes.sort_values('GRANULE_ID')
#
#     novo_df=[]
#     if scenes.empty == False:
#         novo_df = scenes['BASE_URL'].to_list()
#     del scenes, selecao
#
#     cmd = f'/mnt/home/vitorsm/anaconda3/envs/gdaltf/bin/gsutil -q -m cp -r {novo_df[0]} {outfolder}'
#     out1 = os.system(cmd)
#
#     return out1, novo_df[0]
#
# def sort_scenes(scenes):
#     dates = []
#     for scene in scenes:
#         datei = scene.split('T')[0]
#         dates.append((datei, scene))
#     dates.sort()
#
#     scenes = []
#     for i, scene in dates:
#         scenes.append(scene)
#     return scenes
#
# def gdal_to_arraySingleBand(img):
#     '''
#     img: Path to image. This image will be classified. The array is between 0 a 1.
#     :return: array with 3D
#     '''
#     import gdal
#     import numpy as np
#     # read image
#     ds = gdal.Open(img)
#     image = ds.ReadAsArray().astype(np.uint16)  # uint8
#     return np.array(image), ds
#
# def create_temp_folder(temp_repo, k):
#     import string, random, os
#     controli=0
#     while controli==0:
#         letters = string.ascii_lowercase + string.digits
#         result_str = ''.join(random.choice(letters) for i in range(k))
#         temp_folder = os.path.join(temp_repo, 'temp_' + result_str)
#         if not os.path.exists(temp_folder):
#             os.makedirs(temp_folder)
#             controli = 1
#     return temp_folder
#
# def get_analysis_cloud(cloud_tif):
#     import numpy as np
#     checkin = False
#     cloud_arr,ds = gdal_to_arraySingleBand(cloud_tif)
#     r,c = cloud_arr.shape
#     total = r*c
#     n = np.sum(np.ravel(cloud_arr[cloud_arr == 1]))
#     prec = n/float(total)*100
#     if prec > 10: # at least 1% of clouds
#         checkin=True
#
#     return checkin
#
# def create_unmapped(cloudprepost_path, input_water, output_unmapped):
#
#     import numpy as np
#     w_array, ds = gdal_to_arraySingleBand(input_water)
#     c_array, ds = gdal_to_arraySingleBand(cloudprepost_path)
#
#     unmapped = (np.logical_or(c_array == 1, w_array == 1))*1
#     export_tif(output_unmapped, unmapped, ds)
#
# def preproc_cloud(cloud_path, output_path):
#     import numpy as np
#     c_array, ds= gdal_to_arraySingleBand(cloud_path)
#
#     mask_shadow = c_array == 3
#     c_array[mask_shadow]=2
#     mask_clear = c_array != 2
#     c_array[mask_clear] = 0
#     mask_cl = c_array == 2
#     c_array[mask_cl] = 1
#     export_tif(output_path, c_array, ds)
#
#     return
#
# def getdates(datein):
#     data_inicial = f'{datein[0:4]}-{datein[4:6]}-{datein[6:8]}'
#     nextday = int(datein[6:8]) + 1
#     data_final = f'{datein[0:4]}-{datein[4:6]}-{str(nextday).zfill(2)}'
#     return data_inicial, data_final
#
# def write_check(checkpoint):
#     with open(checkpoint, "a+") as myfile:
#         myfile.write(str(1))
#
# def run_fmask(safename, temp_safe):
#     import os
#     import sys
#     sys.path.append(r'Z:\projects\code_center\advance\fmask\python-fmask-0.5.8\fmask')
#     sys.path.append(r'Z:\projects\code_center\advance\fmask\python-fmask-0.5.8\fmask\cmdline')
#     import fmask
#     print(f"Starting Fmask with {safename}")
#     path_safedir1 = os.path.join(temp_safe, safename)
#     nm = safename.replace(".SAFE", "")
#     path_cloud = os.path.join(temp_safe, f"cloud_{nm}.tif")
#
#     path_safedir1=r'Y:\datasets\images\sentinel\level_toa\15SXR\20151121\S2A_MSIL1C_20151121T165612_N0204_R026_T15SXR_20151121T165611.SAFE'
#     path_cloud = fr"C:\Users\vsm71\Downloads\cloud.tif"
#     temp_safe = r'C:\Users\vsm71\Downloads\temp'
#     cmd = rf'C:\anaconda3\envs\testenv\python.exe Z:\projects\code_center\advance\fmask\python-fmask-pythonfmask-0.5.8\bin\fmask_sentinel2Stacked.py --shadowbufferdistance 600 --cloudbufferdistance 300 --pixsize 180 -e {temp_safe} -o {path_cloud} --safedir {path_safedir1}'
#     os.system(cmd)
#
# def run_preproc(temp_1):
#     import glob
#     import os
#
#     clouds = glob.glob(os.path.join(temp_1, "cloud*.tif"))
#     cloudstemps = []
#     for y, cloud_path in enumerate(clouds):
#         output_path_cloud = os.path.join(temp_1, f"CLDPRB{str(y)}.tif")
#         preproc_cloud(cloud_path, output_path_cloud)
#         cloudstemps.append(output_path_cloud)
#
#     return cloudstemps
#
# def findandmoveclouds(safedir1, safedir2, temp_SAFES, temp_1):
#     import shutil
#     sfs = [safedir1, safedir2]
#     for i, sf in enumerate(sfs):
#         nm = sf.replace(".SAFE", "")
#         path_cloud = os.path.join(temp_SAFES, f"cloud_{nm}.tif")
#         outi = os.path.join(temp_1, f"cloud_{str(i)}.tif")
#         shutil.copyfile(path_cloud, outi)
#
# def downloadallSAFES(base_in, REGION, MGRS_TILE, scenes):
#     import os, glob
#     import pandas as pd
#     SAFES_dates=[]
#     SAFES=[]
#     temp_repo='/mnt/ufs18/rs-sat.8/vitorsm/Sentinel_burn/dados/temp'
#     temp_SAFES = create_temp_folder(temp_repo, 15)
#
#     # load the index.csv adapted
#     base_cloude = pd.read_pickle('/mnt/ufs18/rs-sat.6/vitorsm/Sentinel_burn/scripts/download/indexL1.pkl')
#
#     for scene in scenes:
#
#         base_dir_data = os.path.join(base_in, REGION, "T" + MGRS_TILE, scene)
#
#         data_tif = glob.glob(os.path.join(base_dir_data, "data", "combine6B", "*.txt"))[0]
#         cloud_tif = glob.glob(os.path.join(base_dir_data, "data", "cloud_*.tif"))[0]
#         checkin = get_analysis_cloud(cloud_tif)
#         checkin = False  # ------------------------------------------Cancelei o FMASK --------------------
#         if checkin == True:
#             try:
#                 print(f"Starting cloud cover with {base_dir_data}")
#                 _, tilei, datepre, datepost = os.path.basename(data_tif).split(".")[0].split("_")
#
#                 dtds=[datepre, datepost]
#                 for dateii in dtds:
#                     data_inicial, data_final = getdates(dateii)
#
#                     if data_inicial in SAFES_dates:
#                         continue
#
#                     SAFES_dates.append(data_inicial)
#                     out1, pathscene = download_images(base_cloude, MGRS_TILE, data_inicial, data_final, temp_SAFES)
#                     safedir1 = os.path.basename(pathscene)
#                     SAFES.append(safedir1)
#
#             except:
#                 print(f"Error on {scene} {MGRS_TILE}")
#
#         else:
#             print("no cloud or I decide to not use fmask")
#
#     del base_cloude
#
#     return temp_SAFES, SAFES_dates, SAFES
#
#
# def run_processing(scene, MGRS_TILE, REGION, temp_SAFES, SAFES_dates, SAFES):
#
#     """Stack the cloud cover of two images
#
#         Parameters
#         ----------
#         scene : str
#             string with date "20191029"
#         MGRS_TILE : str
#             Name of tile "31PBP"
#         REGION : str
#             Name of region of analysis
#
#         Returns
#         -------
#         None
#         """
#
#     import pandas as pd
#     import gdal
#     import os, shutil, glob
#
#     base_in = f'/mnt/research/{SATPATH}/vitorsm/Sentinel_burn/dados/Classified_pairs'
#     temp_repo = f'/mnt/research/{SATPATH}/vitorsm/Sentinel_burn/dados/temp'
#     base_dir_data = os.path.join(base_in, REGION, "T" + MGRS_TILE, scene)
#
#     cloudcheck = os.path.join(base_dir_data, 'Auxiliary', 'fmask*.txt')
#     if len(cloudcheck) == 1:
#         print(f"Skip: {MGRS_TILE} {scene}")
#         return
#
#     data_tif = glob.glob(os.path.join(base_dir_data, "data", "combine6B", "*.txt"))[0]
#     cloud_tif = os.path.join(base_dir_data, "data", "cloud_prepost.tif")
#     checkin = get_analysis_cloud(cloud_tif)
#     checkin=False  # ------------------------------------------Cancelei o FMASK --------------------
#     if checkin:
#         try:
#             print(f"Starting cloud cover with {base_dir_data}")
#             temp_1 = create_temp_folder(temp_repo, 15)
#             _, tilei, datepre, datepost = os.path.basename(data_tif).split(".")[0].split("_")
#
#             data_inicial_pre, data_final = getdates(datepre)
#             idx = SAFES_dates.index(data_inicial_pre)
#             safedir1 = SAFES[idx]
#             data_inicial_post, data_final = getdates(datepost)
#             idx = SAFES_dates.index(data_inicial_post)
#             safedir2 = SAFES[idx]
#
#             findandmoveclouds(safedir1, safedir2, temp_SAFES, temp_1)
#
#             # Run preprocessing to convert it to binary
#             cloudstemps = run_preproc(temp_1)
#
#             cloudpre = cloudstemps[0]; cloudpost = cloudstemps[1]
#             output_temp = os.path.join(temp_1, 'cloud_temp.tif')
#             cloudprepost_path = os.path.join(base_dir_data, 'data', 'cloud_prepost.tif')
#
#             # rescale it
#             rds = gdal.Open(cloudprepost_path)
#             img_width, img_height = rds.RasterXSize, rds.RasterYSize
#             del rds
#
#             cloudpre_array, cloudpost_array, col, row, isect_bb, ds = intersectraster.runintersect2(cloudpre, cloudpost)
#
#             #----------------------------- pre cloud--------------------------------------------
#             output_temp_pre = os.path.join(temp_1, 'cloud_temp_pre.tif')
#             sentinelcloud.QA_mask_arr_single(cloudpre_array, output_temp_pre, ds)
#             output_temp_pre1 = os.path.join(base_dir_data, 'Auxiliary', 'cpre.tif')
#             sentinelcloud.resample(output_temp_pre, img_height, img_width, output_temp_pre1)
#             #----------------------------- post cloud--------------------------------------------
#             output_temp_post = os.path.join(temp_1, 'cloud_temp_post.tif')
#             sentinelcloud.QA_mask_arr_single(cloudpost_array, output_temp_post, ds)
#             output_temp_post1 = os.path.join(base_dir_data, 'Auxiliary', 'cpost.tif')
#             sentinelcloud.resample(output_temp_post, img_height, img_width, output_temp_post1)
#             #------------------------------ combine----------------------------------------------
#             sentinelcloud.QA_mask_arr(cloudpre_array, cloudpost_array, output_temp, ds)
#
#             os.remove(cloudprepost_path)
#
#             sentinelcloud.resample(output_temp, img_height, img_width, cloudprepost_path)
#
#             shutil.rmtree(temp_1)
#
#             checkpoint = os.path.join(base_dir_data, 'Auxiliary', 'fmask_donecloud.txt')
#             write_check(checkpoint)
#
#         except:
#             print(f"Problem in the download or Fmask processing in the {MGRS_TILE} scene: {scene}")
#             shutil.rmtree(temp_1)
#
#             checkpoint = os.path.join(base_dir_data, 'Auxiliary', 'fmask_errorcloud.txt')
#             write_check(checkpoint)
#     else:
#         checkpoint = os.path.join(base_dir_data, 'Auxiliary', 'fmask_noapplied.txt')
#         write_check(checkpoint)
#
#
#     # MERGE THE WATER MASK AND CLOUD FOR UNMAPPED
#     #cloudprepostfmask_path = glob.glob(os.path.join(base_dir_data, "data", "cloud_*.tif"))[0]
#     #input_water = os.path.join(base_dir_data, 'data', 'watermask.tif')
#     #output_unmapped = os.path.join(base_dir_data, 'data', 'unmapped.tif')
#     #create_unmapped(cloudprepostfmask_path, input_water, output_unmapped)
#
#
# def processing_all_parallel(SAFES, temp_SAFES_l, ncores):
#     with multiprocessing.Pool(processes=ncores) as pool:
#         results=pool.starmap_async(run_fmask, zip(SAFES, temp_SAFES_l)).get()
#         print(results)
#
#
# if __name__ == '__main__':
#
#     import os, glob, sys
#     import multiprocessing
#
#     import time, shutil
#     starttime=time.time()
#
#     MGRS_TILE = ''.join(c for c in str(sys.argv[1]) if c.isprintable())
#     REGION = ''.join(c for c in str(sys.argv[2]) if c.isprintable())
#     SATPATH = ''.join(c for c in str(sys.argv[3]) if c.isprintable())
#
#     base_in = f'/mnt/research/{SATPATH}/vitorsm/Sentinel_burn/dados/Classified_pairs'
#     base_dir = os.path.join(base_in, REGION, "T" + MGRS_TILE)
#     scenes = os.listdir(base_dir)
#     scenes = sort_scenes(scenes)
#     MGRS_TILEs = [MGRS_TILE]*len(scenes)
#     REGIONs = [REGION] * len(scenes)
#
#     #Mudei de ideia sobre o FMASK, e coloco isso aqui
#
#     temp_SAFES, SAFES_dates, SAFES = downloadallSAFES(base_in, REGION, MGRS_TILE, scenes)
#
#     temp_SAFES_l = [temp_SAFES] * len(SAFES)
#
#     #parallel processing
#     ncores = 4
#     processing_all_parallel(SAFES, temp_SAFES_l, ncores)
#
#     for scene in scenes:
#         run_processing(scene, MGRS_TILE, REGION, temp_SAFES, SAFES_dates, SAFES)
#
#     shutil.rmtree(temp_SAFES)
#
#     print("finished")
#     print(starttime-time.time())
