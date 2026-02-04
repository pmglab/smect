import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from tqdm import tqdm
import re
import math
import gzip
from statsmodels.stats.multitest import multipletests

def calculate_LDSC_region_power(folder_path,
                   disease_regions=[(0, 99, 20, 29), (0, 99, 60, 69)],  # Modified to multiple disease regions (x_min, x_max, y_min, y_max)
                   inner_regions=[(0, 99, 14, 35), (0, 99, 54, 75)],    # Modified to multiple inner regions (x_min, x_max, y_min, y_max)
                   sig_threshold=0.05, min_detection_ratio=0.8):
    """
    Calculate statistical power for spatial transcriptomics data (based on multiple rectangular disease regions)
    
    Parameters:
    folder_path: Path to folder containing CSV.GZ files
    disease_regions: List of disease regions, each region is (x_min, x_max, y_min, y_max)
    inner_regions: List of inner regions, each region is (x_min, x_max, y_min, y_max)
    sig_threshold: Significance threshold, default 0.05
    min_detection_ratio: Minimum detection ratio for determining slice as high power, default 0.8 (i.e., 80% points detected as significant)
    
    Returns:
    overall_power: Overall power (average of all slice powers)
    avg_fp_rate: Average false positive rate
    slice_power_details: Detailed power information for each slice
    file_count: Number of successfully processed slices
    """
    
    # Custom coordinate parsing function
    def parse_coordinates(spot):
        """Parse spot format, support C0:0 format"""
        try:
            if pd.isna(spot):
                return None, None
            spot_str = str(spot).strip()
            
            # Process "C0:0" format
            if re.match(r'^C\d+:\d+$', spot_str):
                match = re.match(r'C(\d+):(\d+)', spot_str)
                return int(match.group(1)), int(match.group(2))
        except:
            return None, None
        return None, None
    
    # Check if point is in any of the regions
    def point_in_any_region(x, y, regions):
        """Check if point is in any of the given rectangular regions"""
        for region in regions:
            x_min, x_max, y_min, y_max = region
            if (x_min <= x <= x_max) and (y_min <= y <= y_max):
                return True
        return False
    
    # Calculate total points in regions
    def calculate_total_points(regions):
        """Calculate total points in all regions"""
        total_points = 0
        for region in regions:
            x_min, x_max, y_min, y_max = region
            width = x_max - x_min + 1
            height = y_max - y_min + 1
            total_points += width * height
        return total_points
    
    # Calculate total points in disease regions
    total_disease_points = calculate_total_points(disease_regions)
    
    # Initialize data structures
    slice_powers = []  # Store power for each slice
    all_fp_rates = []  # Store false positive rate for each slice
    slice_power_details = []  # Store detailed information for each slice
    error_files = []
    successful_files = 0
    
    # Get all CSV.GZ files to process
    files = [f for f in os.listdir(folder_path) if f.endswith('.csv.gz')]
    total_files = len(files)
    
    print(f"Total slice count: {total_files}")
    
    # Print disease region information
    print(f"Disease regions ({len(disease_regions)}):")
    for i, region in enumerate(disease_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region {i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (Points: {points})")
    
    # Print inner region information
    print(f"Inner regions ({len(inner_regions)}):")
    for i, region in enumerate(inner_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region {i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (Points: {points})")
    
    print(f"Total disease region points: {total_disease_points}")
    print(f"Minimum detection ratio: {min_detection_ratio}")
    
    # Use tqdm to show progress bar
    for file_name in tqdm(files, desc="Processing files"):
        try:
            # Read CSV.GZ file
            file_path = os.path.join(folder_path, file_name)
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f)
            
            # Add multiple testing correction - use Benjamini-Hochberg method
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
            
            # Count significant points in disease regions for current slice
            disease_significant_count = 0
            disease_total_count = 0
            
            # Count false positive points in outer regions for current slice
            outside_significant_count = 0
            outside_total_count = 0
            
            for i, (x, y) in enumerate(coords):
                x, y = int(x), int(y)
                
                # Check if point is in any disease region
                if point_in_any_region(x, y, disease_regions):
                    disease_total_count += 1
                    if p_values[i] < sig_threshold:
                        disease_significant_count += 1
                
                # Check if point is outside all inner regions (false positive region)
                if not point_in_any_region(x, y, inner_regions):
                    outside_total_count += 1
                    if p_values[i] < sig_threshold:
                        outside_significant_count += 1
            
            # Calculate power for current slice
            if disease_total_count > 0:
                detection_ratio = disease_significant_count / total_disease_points
                # If detection ratio exceeds threshold, power is 1, otherwise calculate actual ratio
                slice_power = 1.0 if detection_ratio >= min_detection_ratio else detection_ratio
            else:
                slice_power = 0.0
                detection_ratio = 0.0
            
            # Calculate false positive rate for current slice (significant point ratio in outer regions)
            fp_rate = outside_significant_count / outside_total_count if outside_total_count > 0 else 0
            
            # Record detailed information for current slice
            slice_info = {
                'file_name': file_name,
                'disease_points_total': disease_total_count,
                'disease_points_significant': disease_significant_count,
                'detection_ratio': detection_ratio,
                'slice_power': slice_power,
                'outside_points_total': outside_total_count,
                'outside_points_significant': outside_significant_count,
                'fp_rate': fp_rate,
                'fp_count': outside_significant_count
            }
            slice_power_details.append(slice_info)
            
            slice_powers.append(slice_power)
            all_fp_rates.append(fp_rate)
            successful_files += 1
            
        except Exception as e:
            print(f"Error processing file {file_name}: {str(e)}")
            error_files.append(file_name)
    
    # Calculate overall metrics
    if successful_files == 0:
        raise ValueError("No valid files processed")
    
    # Overall power = average of all slice powers
    overall_power = np.mean(slice_powers) if slice_powers else 0
    avg_fp_rate = np.mean(all_fp_rates) if all_fp_rates else 0
    
    # Output processing summary
    print(f"\nProcessing summary:")
    print(f"Successfully processed: {successful_files}/{total_files} files")
    print(f"Total disease region points: {total_disease_points}")
    print(f"Overall power: {overall_power:.4f}")
    print(f"Average false positive rate: {avg_fp_rate:.4f}")
    
    if error_files:
        print(f"Failed files: {len(error_files)}")
    
    return overall_power, avg_fp_rate, slice_power_details, successful_files

def save_power_stats_to_csv(slice_power_details, overall_power, avg_fp_rate, file_count, output_folder="C:/Users/pmgl903/Desktop/kgg/kggsum3/simulation/result_ldsc_t400/power_stats"):
    """
    Save power statistics to CSV file
    
    Parameters:
    slice_power_details: Detailed power information for each slice
    overall_power: Overall power
    avg_fp_rate: Average false positive rate
    file_count: Number of successfully processed slices
    output_folder: Path to output folder for CSV files
    """
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # 1. Save overall statistics
    overall_stats = {
        'Overall_power': [overall_power],
        'Average_FP_rate': [avg_fp_rate],
        'Total_slices': [file_count],
        'High_power_slices': [len([s for s in slice_power_details if s['slice_power'] >= 0.8])],
        'Medium_power_slices': [len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8])],
        'Low_power_slices': [len([s for s in slice_power_details if s['slice_power'] < 0.5])],
        'High_power_percentage': [len([s for s in slice_power_details if s['slice_power'] >= 0.8]) / len(slice_power_details) * 100],
        'Medium_power_percentage': [len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]) / len(slice_power_details) * 100],
        'Low_power_percentage': [len([s for s in slice_power_details if s['slice_power'] < 0.5]) / len(slice_power_details) * 100]
    }
    
    overall_df = pd.DataFrame(overall_stats)
    overall_csv_path = os.path.join(output_folder, "overall_power_statistics.csv")
    overall_df.to_csv(overall_csv_path, index=False, encoding='utf-8-sig')
    print(f"Overall statistics saved to: {overall_csv_path}")
    
    # 2. Save detailed information for each slice
    detailed_data = []
    for slice_info in slice_power_details:
        detailed_data.append({
            'Filename': slice_info['file_name'],
            'Disease_significant_points': slice_info['disease_points_significant'],
            'Disease_total_points': slice_info['disease_points_total'],
            'Detection_ratio': slice_info['detection_ratio'],
            'Slice_power': slice_info['slice_power'],
            'Outer_region_points': slice_info['outside_points_total'],
            'False_positive_points': slice_info['fp_count'],
            'False_positive_rate': slice_info['fp_rate']
        })
    
    detailed_df = pd.DataFrame(detailed_data)
    detailed_csv_path = os.path.join(output_folder, "detailed_slice_power_statistics.csv")
    detailed_df.to_csv(detailed_csv_path, index=False, encoding='utf-8-sig')
    print(f"Detailed slice information saved to: {detailed_csv_path}")
    
    # 3. Save power distribution statistics
    power_distribution = {
        'Power_range': ['High_power (≥0.8)', 'Medium_power (0.5-0.8)', 'Low_power (<0.5)'],
        'Slice_count': [
            len([s for s in slice_power_details if s['slice_power'] >= 0.8]),
            len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]),
            len([s for s in slice_power_details if s['slice_power'] < 0.5])
        ],
        'Percentage': [
            len([s for s in slice_power_details if s['slice_power'] >= 0.8]) / len(slice_power_details) * 100,
            len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]) / len(slice_power_details) * 100,
            len([s for s in slice_power_details if s['slice_power'] < 0.5]) / len(slice_power_details) * 100
        ]
    }
    
    distribution_df = pd.DataFrame(power_distribution)
    distribution_csv_path = os.path.join(output_folder, "power_distribution.csv")
    distribution_df.to_csv(distribution_csv_path, index=False, encoding='utf-8-sig')
    print(f"Power distribution statistics saved to: {distribution_csv_path}")
    
    # 4. Save false positive rate statistics
    fp_stats = {
        'Average_FP_rate': [avg_fp_rate],
        'Total_FP_points': [sum([s['fp_count'] for s in slice_power_details])],
        'Average_FP_points_per_slice': [sum([s['fp_count'] for s in slice_power_details]) / len(slice_power_details)],
        'Zero_FP_slices': [len([s for s in slice_power_details if s['fp_count'] == 0])],
        'High_FP_slices(>0.1)': [len([s for s in slice_power_details if s['fp_rate'] > 0.1])]
    }
    
    fp_df = pd.DataFrame(fp_stats)
    fp_csv_path = os.path.join(output_folder, "false_positive_statistics.csv")
    fp_df.to_csv(fp_csv_path, index=False, encoding='utf-8-sig')
    print(f"False positive statistics saved to: {fp_csv_path}")
    
    return overall_csv_path, detailed_csv_path, distribution_csv_path, fp_csv_path


def print_detailed_results(slice_power_details, overall_power, avg_fp_rate, file_count):
    """Print detailed results"""
    print(f"\nDetailed results:")
    print(f"Total slices: {file_count}")
    print(f"Overall power: {overall_power:.4f}")
    print(f"Average false positive rate: {avg_fp_rate:.4f}")
    
    # Sort by power
    sorted_slices = sorted(slice_power_details, key=lambda x: x['slice_power'], reverse=True)
    
    print(f"\nSlice power details (sorted by power descending):")
    print("Filename | Disease_sig/total | Detection_ratio | Slice_power | FP_points | FP_rate")
    print("-" * 100)
    
    for slice_info in sorted_slices[:10]:  # Only show first 10 slices
        print(f"{slice_info['file_name'][:40]}... | {slice_info['disease_points_significant']:3d}/{slice_info['disease_points_total']:3d} | {slice_info['detection_ratio']:8.4f} | {slice_info['slice_power']:8.4f} | {slice_info['fp_count']:6d} | {slice_info['fp_rate']:.4f}")
    
    if len(sorted_slices) > 10:
        print(f"... and {len(sorted_slices) - 10} more slices not shown")
    
    # Count high power slices
    high_power_slices = [s for s in slice_power_details if s['slice_power'] >= 0.8]
    medium_power_slices = [s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]
    low_power_slices = [s for s in slice_power_details if s['slice_power'] < 0.5]
    
    print(f"\nPower statistics:")
    print(f"High power slices(≥0.8): {len(high_power_slices)}/{len(slice_power_details)} ({len(high_power_slices)/len(slice_power_details)*100:.1f}%)")
    print(f"Medium power slices(0.5-0.8): {len(medium_power_slices)}/{len(slice_power_details)} ({len(medium_power_slices)/len(slice_power_details)*100:.1f}%)")
    print(f"Low power slices(<0.5): {len(low_power_slices)}/{len(slice_power_details)} ({len(low_power_slices)/len(slice_power_details)*100:.1f}%)")
    
    # False positive point statistics
    total_fp_points = sum([s['fp_count'] for s in slice_power_details])
    avg_fp_points = total_fp_points / len(slice_power_details) if slice_power_details else 0
    print(f"\nFalse positive point statistics:")
    print(f"Total false positive points: {total_fp_points}")
    print(f"Average false positive points per slice: {avg_fp_points:.2f}")


if __name__ == "__main__":
    # Usage example
    overall_power, fp_rate, slice_details, file_count = calculate_power(
        folder_path=r"C:\Users\pmgl903\Desktop\kgg\kggsum3\simulation\result_ldsc_t400",
        disease_regions=[(20, 39, 20, 39)],  # Two disease regions: X-axis 0-99, Y-axis 20-29 and 60-69
        inner_regions=[(14, 45, 14, 45)], # Two inner regions: X-axis 0-99, Y-axis 14-35 and 54-75
        sig_threshold=0.05,
        min_detection_ratio=0.8  # 80% points detected as significant is considered high power
    )
    print_detailed_results(slice_details, overall_power, fp_rate, file_count)
    
    # New: Save power statistics to CSV files
    csv_files = save_power_stats_to_csv(slice_details, overall_power, fp_rate, file_count)
    
    print(f"\nAll statistics successfully saved to CSV files:")
    for file_path in csv_files:
        print(f"  - {file_path}")