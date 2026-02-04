import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from tqdm import tqdm
import re
import math

def calculate_DESE_onediseasespot_power(folder_path, center=(50,50), radius=3, sigma=1, sig_threshold=0.05):
    """
    Calculate statistical power for spatial transcriptomics data
    
    Parameters:
    folder_path: Folder path containing enrichment.xls files
    center: Center coordinates of disease region, default (30,30)
    radius: Radius of disease region, default 3
    sigma: Sigma parameter for Gaussian decay function, controls weight decay speed with distance, default 1
    sig_threshold: Significance threshold, default 0.05
    
    Returns:
    overall_power: Overall power (average of slice powers)
    avg_fp_rate: Average false positive rate
    slide_details: Detailed information list for each slice
    center_detected_count: Number of slices where center point was detected
    """
    # Define distance weight function (Gaussian decay)
    def distance_weight(distance, sigma=sigma):
        """Calculate distance weight based on Gaussian function"""
        return math.exp(-(distance**2)/(2*sigma**2))
    
    # Custom coordinate parsing function
    def parse_coordinates(condition):
        """Parse coordinate format, supports multiple formats (including C0:X)"""
        try:
            if pd.isna(condition):
                return None, None
            condition_str = str(condition).strip()
            
            # Handle "C0:X" special format
            if condition_str.startswith('C0:'):
                parts = condition_str.split(':')
                if len(parts) >= 2 and parts[1].isdigit():
                    return 0, int(parts[1])
            
            # Handle other formats (C34:30, 34_30, etc.)
            if re.match(r'^C\d+:\d+$', condition_str):
                match = re.match(r'C(\d+):(\d+)', condition_str)
                return int(match.group(1)), int(match.group(2))
            if re.match(r'^\d+_\d+$', condition_str):
                parts = condition_str.split('_')
                return int(parts[0]), int(parts[1])
            if ':' in condition_str:
                parts = condition_str.split(':')
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    return int(parts[0]), int(parts[1])
            if '_' in condition_str:
                parts = condition_str.split('_')
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    return int(parts[0]), int(parts[1])
        except:
            return None, None
        return None, None
    
    # Convert center point to numpy array format
    center = np.array(center).reshape(1, 2)
    disease_radius = radius
    
    # Initialize statistics
    total_files = 0
    slide_powers = []  # Store power for each slice
    all_fp_rates = []  # Store false positive rate for each slice
    error_files = []   # Store filenames that encountered errors
    slide_details = []  # Store detailed information for each slice
    center_detected_count = 0  # Count of slices where center point was detected
    
    # Get all files to process
    files = [f for f in os.listdir(folder_path) 
             if f.startswith('gene.spatial.expression.txt.') 
             and f.endswith('.enrichment.xls')]
    
    # Show progress bar using tqdm
    for file_name in tqdm(files, desc="Processing files"):
        try:
            # Read file
            file_path = os.path.join(folder_path, file_name)
            df = pd.read_excel(file_path, dtype={'Condition': str})
            
            # Parse coordinates
            coordinates = df['Condition'].apply(parse_coordinates)
            df['X'] = [c[0] if c[0] is not None else np.nan for c in coordinates]
            df['Y'] = [c[1] if c[1] is not None else np.nan for c in coordinates]
            df.dropna(subset=['X', 'Y'], inplace=True)
            df['X'] = df['X'].astype(int)
            df['Y'] = df['Y'].astype(int)
            
            # Extract coordinates and p-values
            coords = df[['X', 'Y']].values
            p_values = df['Adjusted(p)'].values
            
            # Calculate distances and weights
            distances = cdist(coords, center)[:, 0]
            in_disease = distances <= disease_radius
            weights = np.array([distance_weight(distances[i]) for i in range(len(distances))]) # Calculate weights for all points
            significant = p_values < sig_threshold
            
            # Check if center point (30,30) exists and is significant
            center_point_mask = (df['X'] == center[0,0]) & (df['Y'] == center[0,1])
            slide_power = 0  # Initialize current slice power
            center_detected = False  # Whether center point was detected
            center_significant = False  # Whether center point is significant
            nearest_sig_dist_in_disease = None # Distance to nearest significant point in disease region
            nearest_sig_weight_in_disease = 0 # Weight of nearest significant point in disease region (will be used as slice power)
            
            if np.any(center_point_mask):
                center_detected = True
                center_detected_count += 1
                center_idx = np.where(center_point_mask)[0][0]
                center_significant = significant[center_idx]
                
                if center_significant:
                    # Center point is significant, slice power is 1
                    slide_power = 1.0
                else:
                    # Center point not significant, need to find nearest significant point in disease region
                    # Get significant points in disease region
                    disease_significant_indices = np.where(in_disease & significant)[0]
                    if disease_significant_indices.size > 0:
                        # Find distance and weight of point closest to center
                        min_dist_index_in_disease = disease_significant_indices[np.argmin(distances[disease_significant_indices])]
                        nearest_sig_dist_in_disease = distances[min_dist_index_in_disease]
                        nearest_sig_weight_in_disease = weights[min_dist_index_in_disease]
                        slide_power = nearest_sig_weight_in_disease
                    else:
                        # No significant points in disease region, power is 0
                        slide_power = 0
            else:
                # Center point doesn't exist, need to find nearest significant point in disease region
                center_detected = False
                # Get significant points in disease region
                disease_significant_indices = np.where(in_disease & significant)[0]
                if disease_significant_indices.size > 0:
                    # Find distance and weight of point closest to center
                    min_dist_index_in_disease = disease_significant_indices[np.argmin(distances[disease_significant_indices])]
                    nearest_sig_dist_in_disease = distances[min_dist_index_in_disease]
                    nearest_sig_weight_in_disease = weights[min_dist_index_in_disease]
                    slide_power = nearest_sig_weight_in_disease
                else:
                    # No significant points in disease region, power is 0
                    slide_power = 0
            
            # Record detailed information for current slice
            slide_info = {
                'file_name': file_name,
                'center_detected': center_detected,
                'center_significant': center_significant if center_detected else None,
                'power': slide_power,
                'nearest_sig_dist_in_disease': nearest_sig_dist_in_disease, # New: record distance to nearest significant point
                'nearest_sig_weight_in_disease': nearest_sig_weight_in_disease # New: record weight of nearest significant point
            }
            slide_details.append(slide_info)
            
            # Record current slice power
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
    
    # Overall power = average of slice powers
    overall_power = np.mean(slide_powers) if slide_powers else 0 # Fixed variable name spelling error slide_power -> slide_powers
    avg_fp_rate = np.mean(all_fp_rates) if all_fp_rates else 0
    
    # Output processing summary
    print(f"\nProcessing Summary:")
    print(f"Successfully processed: {total_files}/{len(files)} files")
    print(f"Number of slices with center point detected: {center_detected_count}/{total_files}")
    if error_files:
        print(f"Files with processing errors: {error_files}")
    
    return overall_power, avg_fp_rate, slide_details, center_detected_count


# if __name__ == "__main__":
#     # Usage example
#     overall_power, fp_rate, slide_details, center_detected_count = calculate_power(
#         folder_path=r"C:\Users\pmgl903\Desktop\kgg\kggsum3\simulation\result_dese_h225",
#         center=(50,50),
#         radius=3,
#         sigma=1.1,
#         sig_threshold=0.05
#     )
    
#     print(f"\nOverall Power: {overall_power:.4f} | False Positive Rate: {fp_rate:.4f}")
#     print(f"Number of slices with center point detected: {center_detected_count}")
    
#     # Output detailed information for each slice
#     print(f"\nDetailed information for each slice:")
#     print("Filename | Center Detected | Center Significant | Power | Nearest Sig Distance | Nearest Sig Weight")
#     print("-" * 100)
#     for detail in slide_details:
#         center_detected = "Yes" if detail['center_detected'] else "No"
#         center_sig = "Yes" if detail['center_significant'] else ("No" if detail['center_significant'] is not None else "N/A")
#         nearest_dist = f"{detail['nearest_sig_dist_in_disease']:.4f}" if detail['nearest_sig_dist_in_disease'] is not None else "N/A"
#         nearest_weight = f"{detail['nearest_sig_weight_in_disease']:.4f}" if detail['nearest_sig_weight_in_disease'] is not None else "N/A"
#         print(f"{detail['file_name']} | {center_detected} | {center_sig} | {detail['power']:.4f} | {nearest_dist} | {nearest_weight}")