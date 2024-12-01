

def extract_landsat(input_gz, output_dir):
    import tarfile, os

    # Create output folder
    os.makedirs(output_dir, exist_ok=True)

    # Extract files
    tar = tarfile.open(input_gz, "r")  # "LC81910182016153-SC20161208043748.tar.gz"
    tar.extractall(output_dir)
    tar.close()