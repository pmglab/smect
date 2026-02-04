import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from tqdm import tqdm
import re
import math
import gzip
from statsmodels.stats.multitest import multipletests

def calculate_LDSC_onediseasespot_power(folder_path, center=(30,30), radius=3, sigma=1, sig_threshold=0.05):
    """
    Calculate statistical power for spatial transcriptomics data

    Parameters:
    folder_path: Path to folder containing CSV.GZ files
    center: Center coordinates of disease region, default (30,30)
    radius: Radius of disease region, default 3
    sigma: Sigma parameter for Gaussian decay function, controls weight decay rate with distance, default 1
    sig_threshold: Significance threshold, default 0.05

    Returns:
    overall_power: Overall power (average of individual slide powers)
    avg_fp_rate: Average false positive rate
    slide_details: List of detailed information for each slide
    center_detected_count: Number of slides where center point was detected
    """
    # Define distance weight function (Gaussian decay)
    def distance_weight(distance, sigma=sigma):
        """Calculate distance weight based on Gaussian function"""
        return math.exp(-(distance**2)/(2*sigma**2))
    
    # Custom coordinate parsing function
    def parse_coordinates(spot):
        """Parse spot format, supports C0:0 format"""
        try:
            if pd.isna(spot):
                return None, None
            spot_str = str(spot).strip()
            
            # Handle "C0:0" format
            if re.match(r'^C\d+:\d+$', spot_str):
                match = re.match(r'C(\d+):(\d+)', spot_str)
                return int(match.group(1)), int(match.group(2))
        except:
            return None, None
        return None, None
    
    # Convert center point to numpy array format
    center = np.array(center).reshape(1, 2)
    disease_radius = radius
    
    # Initialize statistics
    total_files = 0
    slide_powers = []  # Store power for each slide
    all_fp_rates = []  # Store false positive rate for each slide
    error_files = []   # Store filenames that encountered errors during processing
    slide_details = []  # Store detailed information for each slide
    center_detected_count = 0  # Count of slides where center point was detected
    
    # Get all CSV.GZ files to process
    files = [f for f in os.listdir(folder_path) if f.endswith('.csv.gz')]
    
    # Show progress bar using tqdm
    for file_name in tqdm(files, desc="Processing files"):
        try:
            # Read CSV.GZ file
            file_path = os.path.join(folder_path, file_name)
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f)
            
            # Add multiple testing correction - using Benjamini-Hochberg method
            if 'p' in df.columns:
                # Perform FDR correction
                pvals = df['p'].values
                reject, corrected_pvals, _, _ = multipletests(pvals, alpha=sig_threshold, method='fdr_bh')
                df['Adjusted(p)'] = corrected_pvals
            else:
                raise ValueError(f"File {file_name} missing p-value column")
            
            # Parse coordinates
            coordinates = df['spot'].apply(parse_coordinates)
            df['X'] = [c[0] if c[0] is not None else np.nan for c in coordinates]
            df['Y'] = [c[1] if c[1] is not None else np.nan for c in coordinates]
            df.dropna(subset=['X', 'Y'], inplace=True)
            df['X'] = df['X'].astype(int)
            df['Y'] = df['Y'].astype(int)
            
            # Extract coordinates and adjusted p-values
            coords = df[['X', 'Y']].values
            p_values = df['Adjusted(p)'].values
            
            # Calculate distances and weights
            distances = cdist(coords, center)[:, 0]
            in_disease = distances <= disease_radius
            weights = np.array([distance_weight(distances[i]) for i in range(len(distances))])
            significant = p_values < sig_threshold
            
            # Check if center point (30,30) exists and is significant
            center_point_mask = (df['X'] == center[0,0]) & (df['Y'] == center[0,1])
            slide_power = 0  # Initialize power for current slide
            center_detected = False  # Whether center point was detected
            center_significant = False  # Whether center point is significant
            nearest_sig_dist_in_disease = None # Distance to nearest significant point in disease region
            nearest_sig_weight_in_disease = 0 # Weight of nearest significant point in disease region
            
            if np.any(center_point_mask):
                center_detected = True
                center_detected_count += 1
                center_idx = np.where(center_point_mask)[0][0]
                center_significant = significant[center_idx]
                
                if center_significant:
                    # Center point is significant, slide power is 1
                    slide_power = 1.0
                else:
                    # Center point not significant, find nearest significant point in disease region
                    disease_significant_indices = np.where(in_disease & significant)[0]
                    if disease_significant_indices.size > 0:
                        min_dist_index_in_disease = disease_significant_indices[np.argmin(distances[disease_significant_indices])]
                        nearest_sig_dist_in_disease = distances[min_dist_index_in_disease]
                        nearest_sig_weight_in_disease = weights[min_dist_index_in_disease]
                        slide_power = nearest_sig_weight_in_disease
                    else:
                        # No significant points in disease region, power is 0
                        slide_power = 0
            else:
                # Center point doesn't exist, find nearest significant point in disease region
                center_detected = False
                disease_significant_indices = np.where(in_disease & significant)[0]
                if disease_significant_indices.size > 0:
                    min_dist_index_in_disease = disease_significant_indices[np.argmin(distances[disease_significant_indices])]
                    nearest_sig_dist_in_disease = distances[min_dist_index_in_disease]
                    nearest_sig_weight_in_disease = weights[min_dist_index_in_disease]
                    slide_power = nearest_sig_weight_in_disease
                else:
                    # No significant points in disease region, power is 0
                    slide_power = 0
            
            # Record detailed information for current slide
            slide_info = {
                'file_name': file_name,
                'center_detected': center_detected,
                'center_significant': center_significant if center_detected else None,
                'power': slide_power,
                'nearest_sig_dist_in_disease': nearest_sig_dist_in_disease,
                'nearest_sig_weight_in_disease': nearest_sig_weight_in_disease
            }
            slide_details.append(slide_info)
            
            # Record power for current slide
            slide_powers.append(slide_power)
            
            # Calculate false positive rate (proportion of significant points outside disease region)
            outside_disease = ~in_disease
            if np.any(outside_disease):
                fp_rate = np.sum(significant & outside_disease) / np.sum(outside_disease)
                all_fp_rates.append(fp_rate)
            
            total_files += 1
            
        except Exception as e:
            print(f"Error processing file {file_name}: {str(e)}")
            error_files.append(file_name)
    
    # Calculate final metrics
    if total_files == 0:
        raise ValueError("No valid files processed")
    
    # Overall power = average of individual slide powers
    overall_power = np.mean(slide_powers) if slide_powers else 0
    avg_fp_rate = np.mean(all_fp_rates) if all_fp_rates else 0
    
    # Output processing summary
    print(f"\nProcessing Summary:")
    print(f"Successfully processed: {total_files}/{len(files)} files")
    print(f"Slides with center point detected: {center_detected_count}/{total_files}")
    if error_files:
        print(f"Files with processing errors: {error_files}")
    
    return overall_power, avg_fp_rate, slide_details, center_detected_count


# if __name__ == "__main__":
#     # Usage example
#     overall_power, fp_rate, slide_details, center_detected_count = calculate_power(
#         folder_path=r"C:\Users\pmgl903\Desktop\kgg\kggsum3\simulation\result_ldsc_eqtl_h50",
#         center=(50, 50),
#         radius=3,
#         sigma=1.1,
#         sig_threshold=0.05
#     )
    
#     print(f"\nOverall power: {overall_power:.4f} | False positive rate: {fp_rate:.4f}")
#     print(f"Slides with center point detected: {center_detected_count}")
    
#     # Output detailed information for each slide
#     print(f"\nDetailed information for each slide:")
#     print("File Name | Center Detected | Center Significant | Power | Nearest Sig Distance | Nearest Sig Weight")
#     print("-" * 100)
#     for detail in slide_details:
#         center_detected = "Yes" if detail['center_detected'] else "No"
#         center_sig = "Yes" if detail['center_significant'] else ("No" if detail['center_significant'] is not None else "N/A")
#         nearest_dist = f"{detail['nearest_sig_dist_in_disease']:.4f}" if detail['nearest_sig_dist_in_disease'] is not None else "N/A"
#         nearest_weight = f"{detail['nearest_sig_weight_in_disease']:.4f}" if detail['nearest_sig_weight_in_disease'] is not None else "N/A"
#         print(f"{detail['file_name']} | {center_detected} | {center_sig} | {detail['power']:.4f} | {nearest_dist} | {nearest_weight}")