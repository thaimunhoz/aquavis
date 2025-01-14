import geopandas as gpd
import pandas as pd

# Load the Sentinel and Landsat shapefiles
sentinel_tiles = gpd.read_file(r"C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\MGRS_tiles.shp")
landsat_tiles = gpd.read_file(r"C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\ms_landsat.shp")

# Ensure both shapefiles use the same coordinate reference system (CRS)
if sentinel_tiles.crs != landsat_tiles.crs:
    landsat_tiles = landsat_tiles.to_crs(sentinel_tiles.crs)

# Perform spatial intersection
result = []

for _, sentinel_row in sentinel_tiles.iterrows():
    intersecting_landsat = landsat_tiles[landsat_tiles.geometry.intersects(sentinel_row.geometry)]
    landsat_names = intersecting_landsat['Name'].tolist()
    result.append({
        "sentinel_tile": sentinel_row['Name'],
        "sentinel_epsg": sentinel_row['ESPG'],
        "landsat_tiles": landsat_names
    })

# Convert the result to a DataFrame
result_df = pd.DataFrame(result)

# Save the result to a CSV file
result_df.to_csv(r"C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv", index=False)

print("CSV file created successfully!")
