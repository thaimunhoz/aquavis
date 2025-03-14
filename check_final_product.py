import os
import pandas as pd

# Function to check if all required bands are present
def check_bands_exist(hlswater_folder, tile, date):

    required_bands = ['B02', 'B03', 'B04', 'B05']
    missing_bands = []

    # Check for each band
    for band in required_bands:
        band_file = f"HLS_{tile}_{date}_{band}.tif"
        band_path = os.path.join(hlswater_folder, band_file)
        if not os.path.exists(band_path):
            missing_bands.append(band_file)

    return missing_bands


# Function to process each line in the input txt file
def process_files(input_file, output_list):

    tiles_paths = pd.read_csv(input_file)

    for i in range(len(tiles_paths)):

        input_path = tiles_paths['output_path'][i]

        # Extract tile and date information from the input file path
        parts = input_path.split('/')
        tile = parts[-3]  # Assuming tile info is the 3rd to last element
        date = parts[-2]  # Assuming date info is the 2nd to last element

        # Construct the corresponding output path in the 'hlswater' folder
        hlswater_folder = os.path.join("/ddnlus/scratch/r3693/hls_water/HLS_DATASET/sentinel/hlswater", tile)
        output_file_name = f"HLS_{tile}_{date}_S30_v1.0"
        output_path = os.path.join(hlswater_folder, output_file_name)

        # Check if the folder exists
        if not os.path.isdir(hlswater_folder):
            output_list.append(input_path)  # Add input path if folder doesn't exist
        else:
            # Check if the required bands exist
            missing_bands = check_bands_exist(hlswater_folder, tile, date)
            if missing_bands:
                output_list.append(input_path)  # Add input path if bands are missing

# Main function to run the script
def main():
    input_file = r"/ddnlus/r3693/hls_water/scripts/hls_water/src/satwater/auxfiles/tiles/paths_sentinel_toa_35TNH.txt"  # Path to the .txt file with input file paths
    output_list = []  # List to store paths with missing folders or bands

    # Process the files
    process_files(input_file, output_list)

    # Output the results to a file
    with open("missing_files.txt", 'w') as output_file:
        for item in output_list:
            output_file.write(f"{item}\n")

    print(f"Process complete. Missing files listed in 'missing_files.txt'.")


if __name__ == "__main__":
    main()
