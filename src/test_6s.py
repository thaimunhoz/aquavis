import os

os.environ["GDAL_DATA"] = "/ddnlus/r3693/.conda/envs/hls_env/share/gdal"
import rasterio

print(rasterio.env.getenv()["GDAL_DATA"])
