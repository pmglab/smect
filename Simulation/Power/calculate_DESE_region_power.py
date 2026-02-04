import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from tqdm import tqdm
import re
import math
from collections import defaultdict

def calculate_DESE_region_power(folder_path,
                   disease_regions=[(20, 29, 0, 99), (60, 69, 0, 99)],  # Modified to multiple disease regions
                   inner_regions=[(14, 35, 0, 99), (54, 75, 0, 99)],    # Modified to multiple inner regions
                   power_threshold=0.8, sig_threshold=0.05):
    """
    Calculate statistical power for spatial transcriptomics data (based on multiple rectangular disease regions)
    
    Parameters:
    folder_path: Folder path containing enrichment.xls files
    disease_regions: List of disease regions, each region is (x_min, x_max, y_min, y_max)
    inner_regions: List of inner regions, each region is (x_min, x_max, y_min, y_max)
    power_threshold: Threshold ratio to determine high power slices, default 0.8 (80% points significant)
    sig_threshold: Significance threshold, default 0.05
    
    Returns:
    overall_power: Overall power (average of all slice powers)
    avg_fp_rate: Average false positive rate
    slice_power_details: Detailed power information for each slice
    total_files: Number of successfully processed slices
    """
    
    # Custom coordinate parsing function
    def parse_coordinates(condition):
        """Parse coordinate format, supporting multiple formats (including C0:X)"""
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
    
    # Calculate total points in disease regions (sum of all regions) [6,7](@ref)
    def calculate_total_disease_points(regions):
        """Calculate total points in all disease regions"""
        total_points = 0
        for region in regions:
            x_min, x_max, y_min, y_max = region
            width = x_max - x_min + 1
            height = y_max - y_min + 1
            total_points += width * height
        return total_points
    
    # Check if point is in any of the given regions [6,7](@ref)
    def is_point_in_any_region(x, y, regions):
        """Check if point is in any of the given rectangular regions"""
        for region in regions:
            x_min, x_max, y_min, y_max = region
            if (x_min <= x <= x_max) and (y_min <= y <= y_max):
                return True
        return False
    
    # Initialize statistics
    total_files = 0
    error_files = []   # Store files that encountered processing errors
    
    # Get all files to process
    files = [f for f in os.listdir(folder_path) 
             if f.startswith('gene.spatial.expression.txt.') 
             and f.endswith('.enrichment.xls')]
    
    # Calculate total points in disease regions
    total_disease_points = calculate_total_disease_points(disease_regions)
    
    # Store information for each slice
    slice_power_details = []
    
    # Use tqdm to show progress bar
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
            
            # Count significant points in disease regions for current slice
            disease_significant_count = 0
            disease_total_count = 0
            
            # Count false positive points outside inner regions for current slice
            outside_significant_count = 0
            outside_total_count = 0
            
            for _, row in df.iterrows():
                x, y = int(row['X']), int(row['Y'])
                
                # Check if point is in any disease region [6,7](@ref)
                if is_point_in_any_region(x, y, disease_regions):
                    disease_total_count += 1
                    if row['Adjusted(p)'] < sig_threshold:
                        disease_significant_count += 1
                
                # Check if point is outside all inner regions (false positive region) [6,7](@ref)
                if not is_point_in_any_region(x, y, inner_regions):
                    outside_total_count += 1
                    if row['Adjusted(p)'] < sig_threshold:
                        outside_significant_count += 1
            
            # Calculate power for current slice
            if disease_total_count > 0:
                # Power calculation: ratio of significant points, power is 1 if above threshold
                detection_ratio = disease_significant_count / total_disease_points
                slice_power = 1.0 if detection_ratio >= power_threshold else detection_ratio
            else:
                slice_power = 0.0
                detection_ratio = 0.0
            
            # Calculate false positive rate for current slice
            fp_rate = outside_significant_count / outside_total_count if outside_total_count > 0 else 0
            
            # Store current slice information
            slice_info = {
                'file_name': file_name,
                'disease_points_detected': disease_total_count,
                'disease_points_significant': disease_significant_count,
                'disease_points_total': total_disease_points,
                'detection_ratio': detection_ratio,
                'slice_power': slice_power,
                'outside_points_total': outside_total_count,
                'outside_points_significant': outside_significant_count,
                'fp_rate': fp_rate,
                'fp_count': outside_significant_count  # Number of false positive points
            }
            
            slice_power_details.append(slice_info)
            total_files += 1
            
        except Exception as e:
            print(f"Error processing file {file_name}: {str(e)}")
            error_files.append(file_name)
    
    if total_files == 0:
        raise ValueError("No valid files processed")
    
    # Calculate overall power (average of all slice powers)
    overall_power = np.mean([slice_info['slice_power'] for slice_info in slice_power_details])
    
    # Calculate average false positive rate
    avg_fp_rate = np.mean([slice_info['fp_rate'] for slice_info in slice_power_details])
    
    # Print processing summary
    print(f"\nProcessing Summary:")
    print(f"Successfully processed: {total_files}/{len(files)} files")
    
    # Print disease region information [6,7](@ref)
    print(f"Disease Regions ({len(disease_regions)}):")
    for i, region in enumerate(disease_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region{i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (Points: {points})")
    
    # Print inner region information [6,7](@ref)
    print(f"Inner Regions ({len(inner_regions)}):")
    for i, region in enumerate(inner_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region{i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (Points: {points})")
    
    print(f"Total disease region points: {total_disease_points}")
    print(f"Power threshold: {power_threshold} ({power_threshold*100}% points significant)")
    if error_files:
        print(f"Files with processing errors: {error_files}")
    
    return overall_power, avg_fp_rate, slice_power_details, total_files

def save_power_stats_to_csv(slice_power_details, overall_power, avg_fp_rate, file_count, output_folder="C:/Users/pmgl903/Desktop/kgg/kggsum3/simulation/result_DSEE_3568/power_stats/power_stats"):
    """
    Save power statistics to CSV files [1,2](@ref)
    
    Parameters:
    slice_power_details: Detailed power information for each slice
    overall_power: Overall power
    avg_fp_rate: Average false positive rate
    file_count: Number of successfully processed slices
    output_folder: Folder path for output CSV files
    """
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # 1. Save overall statistics [1](@ref)
    overall_stats = {
        'Overall Power': [overall_power],
        'Average False Positive Rate': [avg_fp_rate],
        'Total Slices': [file_count],
        'High Power Slices': [len([s for s in slice_power_details if s['slice_power'] >= 0.8])],
        'Medium Power Slices': [len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8])],
        'Low Power Slices': [len([s for s in slice_power_details if s['slice_power'] < 0.5])],
        'High Power Slice Percentage': [len([s for s in slice_power_details if s['slice_power'] >= 0.8]) / len(slice_power_details) * 100],
        'Medium Power Slice Percentage': [len([s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]) / len(slice_power_details) * 100],
        'Low Power Slice Percentage': [len([s for s in slice_power_details if s['slice_power'] < 0.5]) / len(slice_power_details) * 100]
    }
    
    overall_df = pd.DataFrame(overall_stats)
    overall_csv_path = os.path.join(output_folder, "overall_power_statistics.csv")
    overall_df.to_csv(overall_csv_path, index=False, encoding='utf-8-sig')
    print(f"Overall statistics saved to: {overall_csv_path}")
    
    # 2. Save detailed information for each slice [1,3](@ref)
    detailed_data = []
    for slice_info in slice_power_details:
        detailed_data.append({
            'Filename': slice_info['file_name'],
            'Significant Disease Points': slice_info['disease_points_significant'],
            'Total Disease Points': slice_info['disease_points_total'],
            'Detection Ratio': slice_info['detection_ratio'],
            'Slice Power': slice_info['slice_power'],
            'Total Points Outside': slice_info['outside_points_total'],
            'False Positive Points': slice_info['fp_count'],
            'False Positive Rate': slice_info['fp_rate']
        })
    
    detailed_df = pd.DataFrame(detailed_data)
    detailed_csv_path = os.path.join(output_folder, "detailed_slice_power_statistics.csv")
    detailed_df.to_csv(detailed_csv_path, index=False, encoding='utf-8-sig')
    print(f"Detailed slice information saved to: {detailed_csv_path}")
    
    # 3. Save power distribution statistics [1](@ref)
    power_distribution = {
        'Power Range': ['High Power (≥0.8)', 'Medium Power (0.5-0.8)', 'Low Power (<0.5)'],
        'Number of Slices': [
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
    
    # 4. Save false positive statistics [1](@ref)
    fp_stats = {
        'Average False Positive Rate': [avg_fp_rate],
        'Total False Positive Points': [sum([s['fp_count'] for s in slice_power_details])],
        'Average False Positive Points per Slice': [sum([s['fp_count'] for s in slice_power_details]) / len(slice_power_details)],
        'Slices with Zero False Positives': [len([s for s in slice_power_details if s['fp_count'] == 0])],
        'High False Positive Slices(>0.1)': [len([s for s in slice_power_details if s['fp_rate'] > 0.1])]
    }
    
    fp_df = pd.DataFrame(fp_stats)
    fp_csv_path = os.path.join(output_folder, "false_positive_statistics.csv")
    fp_df.to_csv(fp_csv_path, index=False, encoding='utf-8-sig')
    print(f"False positive statistics saved to: {fp_csv_path}")
    
    return overall_csv_path, detailed_csv_path, distribution_csv_path, fp_csv_path

def print_detailed_results(slice_power_details, overall_power, avg_fp_rate, file_count):
    """Print detailed results"""
    print(f"\nDetailed Results:")
    print(f"Total slices: {file_count}")
    print(f"Overall power: {overall_power:.4f}")
    print(f"Average false positive rate: {avg_fp_rate:.4f}")
    
    # Sort by power
    sorted_slices = sorted(slice_power_details, key=lambda x: x['slice_power'], reverse=True)
    
    print(f"\nSlice power details (sorted by power, descending):")
    print("Filename | Significant/Total Disease Points | Detection Ratio | Slice Power | FP Points | FP Rate")
    print("-" * 100)
    
    for slice_info in sorted_slices[:10]:  # Only show first 10 slices
        print(f"{slice_info['file_name'][:40]}... | {slice_info['disease_points_significant']:3d}/{slice_info['disease_points_total']:3d} | {slice_info['detection_ratio']:8.4f} | {slice_info['slice_power']:8.4f} | {slice_info['fp_count']:6d} | {slice_info['fp_rate']:.4f}")
    
    if len(sorted_slices) > 10:
        print(f"... and {len(sorted_slices) - 10} more slices not shown")
    
    # Count high power slices
    high_power_slices = [s for s in slice_power_details if s['slice_power'] >= 0.8]
    medium_power_slices = [s for s in slice_power_details if 0.5 <= s['slice_power'] < 0.8]
    low_power_slices = [s for s in slice_power_details if s['slice_power'] < 0.5]
    
    print(f"\nPower Statistics:")
    print(f"High power slices(≥0.8): {len(high_power_slices)}/{len(slice_power_details)} ({len(high_power_slices)/len(slice_power_details)*100:.1f}%)")
    print(f"Medium power slices(0.5-0.8): {len(medium_power_slices)}/{len(slice_power_details)} ({len(medium_power_slices)/len(slice_power_details)*100:.1f}%)")
    print(f"Low power slices(<0.5): {len(low_power_slices)}/{len(slice_power_details)} ({len(low_power_slices)/len(slice_power_details)*100:.1f}%)")
    
    # False positive point statistics
    total_fp_points = sum([s['fp_count'] for s in slice_power_details])
    avg_fp_points = total_fp_points / len(slice_power_details) if slice_power_details else 0
    print(f"\nFalse Positive Point Statistics:")
    print(f"Total false positive points: {total_fp_points}")
    print(f"Average false positive points per slice: {avg_fp_points:.2f}")

# if __name__ == "__main__":
#     # Usage example
#     overall_power, fp_rate, slice_details, file_count = calculate_power(
#         folder_path=r"C:\Users\pmgl903\Desktop\kgg\kggsum3\simulation\result_DSEE_3568",
#         disease_regions=[(30, 39, 30, 39), (50,59,50,59), (60,69,60,69),(80,89,80,89)], # Two disease regions
#         inner_regions=[(24, 45, 24, 45), (44,65,44,65), (54,75,54,75),(74,95,74,95)],
#         power_threshold=0.8,  # Slice power is 1 when 80% points are significant
#         sig_threshold=0.05
#     )
#     print_detailed_results(slice_details, overall_power, fp_rate, file_count)
    
#     # New: Save power statistics to CSV files [1,2](@ref)
#     csv_files = save_power_stats_to_csv(slice_details, overall_power, fp_rate, file_count)
    
#     print(f"\nAll statistics successfully saved to CSV files:")
#     for file_path in csv_files:
#         print(f"  - {file_path}")