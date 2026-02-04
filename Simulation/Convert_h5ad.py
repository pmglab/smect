import os
import pandas as pd
import numpy as np
import anndata as ad
from glob import glob

def convert_spatial_to_h5ad(input_dir, output_subdir="h5ad_output", region_ranges=None):
    """
    Convert spatial transcriptomics text data to h5ad format, supporting user-specified region ranges.
    
    Parameters:
        input_dir (str): Input directory path
        output_subdir (str): Output subdirectory name, default "h5ad_output"
        region_ranges (list): List of region ranges, each element is a dict containing region name and coordinate ranges
            Example: [{"name": "Layer1", "x_range": (27, 30), "y_range": (27, 30)}]
            If None, use default ranges
    
    Returns:
        str: Output directory path
    """
    # Set default region ranges
    if region_ranges is None:
        region_ranges = [
            {"name": "Layer1", "x_range": (27, 30), "y_range": (27, 30)},
            {"name": "Layer2", "x_range": None, "y_range": None}  # Default region
        ]
    
    # Set output path
    output_dir = os.path.join(input_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all matching files
    file_pattern = os.path.join(input_dir, "gene.spatial.expression.txt.*")
    file_paths = glob(file_pattern)
    
    if not file_paths:
        print(f"No matching files found in {input_dir}")
        return output_dir
    
    # Process each file
    for file_path in file_paths:
        try:
            adata = process_single_file(file_path, region_ranges)
            save_h5ad_file(adata, file_path, output_dir)
            print(f"Successfully converted: {os.path.basename(file_path)}")
            
            # Print region statistics
            print_region_statistics(adata, region_ranges)
            
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
            continue
    
    print(f"Conversion completed! Results saved in: {output_dir}")
    return output_dir


def process_single_file(file_path, region_ranges):
    """
    Process a single file and convert it to an AnnData object.
    
    Parameters:
        file_path (str): Input file path
        region_ranges (list): List of region ranges
    
    Returns:
        ad.AnnData: Converted AnnData object
    """
    # Read data
    df = pd.read_csv(file_path, sep="\t")
    
    # Process column names and gene names
    df, gene_names = extract_gene_names(df)
    
    # Transpose data and set gene names as column names
    df = df.T
    df.columns = gene_names
    
    # Parse coordinates
    coords = parse_coordinates(df.index)
    
    # Create observation names (add C prefix)
    obs_names = create_observation_names(df.index)
    
    # Create AnnData object
    adata = create_anndata_object(df, obs_names, gene_names)
    
    # Add spatial coordinates
    adata.obsm["spatial"] = coords
    adata.layers["count"] = adata.X
    
    # Add annotations based on user-specified regions
    add_custom_annotations(adata, coords, region_ranges)
    
    return adata


def add_custom_annotations(adata, coords, region_ranges):
    """
    Add annotations based on user-specified region ranges.
    
    Parameters:
        adata (ad.AnnData): AnnData object
        coords (np.array): Coordinate array
        region_ranges (list): List of region ranges
    """
    annotations = []
    
    for coord in coords:
        x, y = coord
        assigned_region = None
        
        # Check if the point is in any specified region
        for region in region_ranges:
            if region["name"] == "Layer2" and region["x_range"] is None:
                # Default region, used when not matching any specific region
                continue
                
            x_range = region["x_range"]
            y_range = region["y_range"]
            
            # Check if the point is within the current region
            if (x_range[0] <= x <= x_range[1]) and (y_range[0] <= y <= y_range[1]):
                assigned_region = region["name"]
                break
        
        # If not matching any specific region, assign to default region
        if assigned_region is None:
            # Find the default region (usually the last one)
            for region in region_ranges:
                if region["x_range"] is None or region["name"] == "Layer2":
                    assigned_region = region["name"]
                    break
        
        annotations.append(assigned_region)
    
    adata.obs["annotation"] = pd.Categorical(annotations)


def is_point_in_region(x, y, region):
    """
    Check if a point is within a specified region.
    
    Parameters:
        x (float): X-coordinate of the point
        y (float): Y-coordinate of the point
        region (dict): Region information, containing x_range and y_range
    
    Returns:
        bool: Whether the point is within the region
    """
    if region["x_range"] is None or region["y_range"] is None:
        return False
    
    x_range = region["x_range"]
    y_range = region["y_range"]
    
    return (x_range[0] <= x <= x_range[1]) and (y_range[0] <= y <= y_range[1])


def print_region_statistics(adata, region_ranges):
    """
    Print statistical information for each region.
    
    Parameters:
        adata (ad.AnnData): AnnData object
        region_ranges (list): List of region ranges
    """
    print("Region statistics:")
    for region in region_ranges:
        region_name = region["name"]
        count = (adata.obs["annotation"] == region_name).sum()
        print(f"  {region_name}: {count} points")
        
        # Print region range information
        if region["x_range"] is not None:
            x_range = region["x_range"]
            y_range = region["y_range"]
            print(f"    Range: x={x_range[0]}-{x_range[1]}, y={y_range[0]}-{y_range[1]}")


# The following helper functions remain unchanged (same as before)
def extract_gene_names(df):
    """Extract gene names from the data frame"""
    if df.columns[0] == "gene":
        gene_names = df["gene"].values
        df = df.drop(columns="gene")
    else:
        gene_names = df.iloc[:, 0].values
        df = df.drop(df.columns[0], axis=1)
    
    return df, gene_names


def parse_coordinates(index):
    """Parse coordinate information from index strings"""
    coords = []
    for idx in index:
        coord_parts = idx.split(":")
        if len(coord_parts) >= 2:
            coords.append([int(coord_parts[0]), int(coord_parts[1])])
        else:
            coords.append([0, 0])
    
    return np.array(coords)


def create_observation_names(index):
    """Create observation names (add C prefix)"""
    return ["C" + str(idx) for idx in index]


def create_anndata_object(df, obs_names, gene_names):
    """Create AnnData object"""
    return ad.AnnData(
        X=df.astype(np.float32).values,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=df.columns)
    )


def save_h5ad_file(adata, file_path, output_dir):
    """Save h5ad file"""
    base_name = os.path.basename(file_path)
    output_path = os.path.join(output_dir, base_name + ".h5ad")
    adata.write_h5ad(output_path)


# Usage example
if __name__ == "__main__":
    input_directory = "./Sim_hot30/ExpressionPhenotypeSimulateTask/"
    
    # Example 1: Use default ranges (27<=x<=30 and 27<=y<=30)
    print("=== Using default ranges ===")
    convert_spatial_to_h5ad(input_directory)
    
    # Example 2: Use custom ranges
    print("\n=== Using custom ranges ===")
    custom_ranges = [
        {"name": "Hotspot1", "x_range": (25, 28), "y_range": (25, 28)},
        {"name": "Hotspot2", "x_range": (32, 35), "y_range": (32, 35)},
        {"name": "Background", "x_range": None, "y_range": None}  # Default region
    ]
    convert_spatial_to_h5ad(input_directory, "h5ad_custom", custom_ranges)
    
    # Example 3: Use multiple complex regions
    print("\n=== Using multiple complex regions ===")
    complex_ranges = [
        {"name": "Region_A", "x_range": (20, 25), "y_range": (20, 25)},
        {"name": "Region_B", "x_range": (30, 35), "y_range": (30, 35)},
        {"name": "Region_C", "x_range": (40, 45), "y_range": (40, 45)},
        {"name": "Other", "x_range": None, "y_range": None}
    ]
    convert_spatial_to_h5ad(input_directory, "h5ad_complex", complex_ranges)
