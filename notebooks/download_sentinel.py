# Utilities
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import requests

from sentinelhub import SHConfig, DataCollection, SentinelHubCatalog, SentinelHubRequest, BBox, bbox_to_dimensions, CRS, MimeType, Geometry

# Set up Sentinel Hub configuration
config = SHConfig()
config.sh_client_id = "sh-dc9d4956-9bff-4dd6-93fc-f87de0294de1"
config.sh_client_secret = "oQn19VXxoQXfyymSQvNGxtIAxiX1IHg9"

config.sh_token_url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
config.sh_base_url = "https://sh.dataspace.copernicus.eu"
config.save("cdse")
# Saved config can be later accessed with config = SHConfig("cdse")

#config = SHConfig("cdse")

catalog = SentinelHubCatalog(config=config)

time_interval = '2023-01-01', '2023-12-31'
aoi_bbox = BBox(bbox=[-58.088837,-34.909584,-57.729034,-34.617387], crs=CRS.WGS84)

search_iterator = catalog.search(
    DataCollection.SENTINEL2_L1C,
    bbox=aoi_bbox,
    time=time_interval,
    fields={"include": ["id", "properties.datetime"], "exclude": []},

)

results = list(search_iterator)
print("Total number of results:", len(results))

evalscript_all_bands = """
    //VERSION=3
    function setup() {
        return {
            input: [{
                bands: ["B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B10","B11","B12"],
                units: "DN"
            }],
            output: {
                bands: 13,
                sampleType: "INT16"
            }
        };
    }

    function evaluatePixel(sample) {
        return [sample.B01,
                sample.B02,
                sample.B03,
                sample.B04,
                sample.B05,
                sample.B06,
                sample.B07,
                sample.B08,
                sample.B8A,
                sample.B09,
                sample.B10,
                sample.B11,
                sample.B12];
    }
"""
request_all_bands = SentinelHubRequest(
    data_folder=r"Z:\guser\tml\mypapers\HLS_package_paper\sentinel_toa",
    input_data=[
        SentinelHubRequest.input_data(
            data_collection=DataCollection.SENTINEL2_L1C,
            time_interval=("2023-01-01", "2023-12-31"),
            maxcc=0.3

        )
    ],
    responses=[SentinelHubRequest.output_response("default", MimeType.tif)],
    bbox=aoi_bbox,
    config=config,
)

all_bands_img = request_all_bands.save_data()

download_dir = r"Z:\guser\tml\mypapers\HLS_package_paper\sentinel_toa"
os.makedirs(download_dir, exist_ok=True)

for result in results:

    product_id = result['id']

    if '21HVB' not in product_id:
        continue

    print(f"Processing product ID: {product_id}")

    # Get product info
    product_info = catalog.get_product(product_id, DataCollection.SENTINEL2_L1C)
    download_url = product_info['assets']['data']['href']  # URL for the .SAFE file

    # Download the .SAFE file
    output_path = os.path.join(download_dir, f"{product_id}.SAFE")
    with requests.get(download_url, stream=True) as response:
        if response.status_code == 200:
            with open(output_path, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
            print(f"Downloaded: {output_path}")
        else:
            print(f"Failed to download product {product_id}. HTTP status: {response.status_code}")