import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from tqdm import tqdm
import re
import math
import gzip
from statsmodels.stats.multitest import multipletests

def calculate_scDRS_region_power(folder_path,
                   disease_regions=[(30, 39, 30, 39), (50, 59, 50, 59),(60, 69, 60, 69)],  # Multiple disease regions (x_min, x_max, y_min, y_max)
                   inner_regions=[(24, 45, 24, 45),(44, 65, 44, 65),(54, 75, 54, 75)],    # Multiple inner regions (x_min, x_max, y_min, y_max)
                   power_threshold=0.8, sig_threshold=0.05):
    """
    Calculate statistical power for spatial transcriptomics data (based on multiple rectangular disease regions)
    
    Parameters:
    folder_path: Path to folder containing .score.gz files
    disease_regions: List of disease regions, each as (x_min, x_max, y_min, y_max)
    inner_regions: List of inner regions, each as (x_min, x_max, y_min, y_max)
    power_threshold: Threshold for determining slide as high-power, default 0.8 (80% points significant)
    sig_threshold: Significance threshold, default 0.05
    
    Returns:
    overall_power: Overall power (average of all slide powers)
    avg_fp_rate: Average false positive rate
    slide_details: List of detailed information for each slide
    total_files: Number of successfully processed slides
    """
    
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
    
    # Check if point is in any of the given regions
    def point_in_any_region(x, y, regions):
        """Check if point is within any of the given rectangular regions"""
        for region in regions:
            x_min, x_max, y_min, y_max = region
            if (x_min <= x <= x_max) and (y_min <= y <= y_max):
                return True
        return False
    
    # Calculate total points in all regions
    def calculate_total_points(regions):
        """Calculate total number of points across all regions"""
        total_points = 0
        for region in regions:
            x_min, x_max, y_min, y_max = region
            width = x_max - x_min + 1
            height = y_max - y_min + 1
            total_points += width * height
        return total_points
    
    # Calculate total points in disease regions
    total_disease_points = calculate_total_points(disease_regions)
    
    # Initialize statistics
    total_files = 0
    slide_powers = []  # Store power for each slide
    all_fp_rates = []  # Store false positive rate for each slide
    error_files = []   # Store filenames that encountered errors
    slide_details = []  # Store detailed information for each slide
    
    # Get all .score.gz files to process
    files = [f for f in os.listdir(folder_path) if f.endswith('.score.gz')]
    
    print(f"Total slides: {len(files)}")
    
    # Print disease region information
    print(f"Disease regions ({len(disease_regions)}):")
    for i, region in enumerate(disease_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region {i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (points: {points})")
    
    # Print inner region information
    print(f"Inner regions ({len(inner_regions)}):")
    for i, region in enumerate(inner_regions):
        x_min, x_max, y_min, y_max = region
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        points = width * height
        print(f"  Region {i+1}: X[{x_min}-{x_max}], Y[{y_min}-{y_max}] (points: {points})")
    
    print(f"Total disease points: {total_disease_points}")
    print(f"Power threshold: {power_threshold} ({power_threshold*100}% points significant)")
    
    # Use tqdm to show progress bar
    for file_name in tqdm(files, desc="Processing files"):
        try:
            # Read .score.gz file
            file_path = os.path.join(folder_path, file_name)
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t')
            
            # Add multiple testing correction - using Benjamini-Hochberg method (using pval column)
            if 'pval' in df.columns:
                # Perform FDR correction
                pvals = df['pval'].values
                reject, corrected_pvals, _, _ = multipletests(pvals, alpha=sig_threshold, method='fdr_bh')
                df['Adjusted(p)'] = corrected_pvals
            else:
                raise ValueError(f"File {file_name} missing pval column")
            
            # Parse coordinates
            coordinates = df['Unnamed: 0'].apply(parse_coordinates)
            df['X'] = [c[0] if c[0] is not None else np.nan for c in coordinates]
            df['Y'] = [c[1] if c[1] is not None else np.nan for c in coordinates]
            df.dropna(subset=['X', 'Y'], inplace=True)
            df['X'] = df['X'].astype(int)
            df['Y'] = df['Y'].astype(int)
            
            # Extract coordinates and corrected p-values
            coords = df[['X', 'Y']].values
            p_values = df['Adjusted(p)'].values
            
            # Count significant points in disease regions for current slide
            disease_significant_count = 0
            disease_detected_count = 0
            
            # Count false positive points outside inner regions for current slide
            outside_significant_count = 0
            outside_total_count = 0
            
            for i, (x, y) in enumerate(coords):
                x, y = int(x), int(y)
                
                # Check if point is in any disease region
                if point_in_any_region(x, y, disease_regions):
                    disease_detected_count += 1
                    if p_values[i] < sig_threshold:
                        disease_significant_count += 1
                
                # Check if point is outside all inner regions (false positive region)
                if not point_in_any_region(x, y, inner_regions):
                    outside_total_count += 1
                    if p_values[i] < sig_threshold:
                        outside_significant_count += 1
            
            # Calculate power for current slide
            if total_disease_points > 0:
                detection_ratio = disease_significant_count / total_disease_points
                # If detection ratio exceeds threshold, power is 1, otherwise use actual ratio
                slide_power = 1.0 if detection_ratio >= power_threshold else detection_ratio
            else:
                slide_power = 0.0
                detection_ratio = 0.0
            
            # Calculate false positive rate for current slide
            fp_rate = outside_significant_count / outside_total_count if outside_total_count > 0 else 0
            
            # Record detailed information for current slide
            slide_info = {
                'file_name': file_name,
                'disease_points_detected': disease_detected_count,
                'disease_points_significant': disease_significant_count,
                'disease_points_total': total_disease_points,
                'detection_ratio': detection_ratio,
                'power': slide_power,
                'outside_points_total': outside_total_count,
                'outside_points_significant': outside_significant_count,
                'fp_rate': fp_rate,
                'fp_count': outside_significant_count  # Number of false positive points
            }
            
            slide_details.append(slide_info)
            slide_powers.append(slide_power)
            all_fp_rates.append(fp_rate)
            total_files += 1
            
        except Exception as e:
            print(f"Error processing file {file_name}: {str(e)}")
            error_files.append(file_name)
    
    # Calculate final metrics
    if total_files == 0:
        raise ValueError("No valid files processed")
    
    # Overall power = average of all slide powers
    overall_power = np.mean(slide_powers) if slide_powers else 0
    avg_fp_rate = np.mean(all_fp_rates) if all_fp_rates else 0
    
    # Output processing summary
    print(f"\nProcessing summary:")
    print(f"Successfully processed: {total_files}/{len(files)} files")
    print(f"Total disease points: {total_disease_points}")
    print(f"Overall power: {overall_power:.4f}")
    print(f"Average false positive rate: {avg_fp_rate:.4f}")
    if error_files:
        print(f"Failed files: {len(error_files)}")
    
    return overall_power, avg_fp_rate, slide_details, total_files


def print_detailed_results(slide_details, overall_power, avg_fp_rate, file_count):
    """Print detailed results"""
    print(f"\nDetailed results:")
    print(f"Total slides: {file_count}")
    print(f"Overall power: {overall_power:.4f}")
    print(f"Average false positive rate: {avg_fp_rate:.4f}")
    
    # Sort by power
    sorted_slides = sorted(slide_details, key=lambda x: x['power'], reverse=True)
    
    print(f"\nSlide power details (sorted by power descending, showing top 10):")
    print("File name | Disease points significant/total | Detection ratio | Slide power | FP points | FP rate")
    print("-" * 100)
    
    for slide_info in sorted_slides[:100]:
        print(f"{slide_info['file_name'][:40]}... | {slide_info['disease_points_significant']:3d}/{slide_info['disease_points_total']:3d} | {slide_info['detection_ratio']:8.4f} | {slide_info['power']:8.4f} | {slide_info['fp_count']:6d} | {slide_info['fp_rate']:.4f}")
    
    if len(sorted_slides) > 10:
        print(f"... {len(sorted_slides) - 10} more slides not shown")
    
    # Count high-power slides
    high_power_slides = [s for s in slide_details if s['power'] >= 0.8]
    medium_power_slides = [s for s in slide_details if 0.5 <= s['power'] < 0.8]
    low_power_slides = [s for s in slide_details if s['power'] < 0.5]
    
    print(f"\nPower statistics:")
    print(f"High-power slides(≥0.8): {len(high_power_slides)}/{len(slide_details)} ({len(high_power_slides)/len(slide_details)*100:.1f}%)")
    print(f"Medium-power slides(0.5-0.8): {len(medium_power_slides)}/{len(slide_details)} ({len(medium_power_slides)/len(slide_details)*100:.1f}%)")
    print(f"Low-power slides(<0.5): {len(low_power_slides)}/{len(slide_details)} ({len(low_power_slides)/len(slide_details)*100:.1f}%)")
    
    # False positive point statistics
    total_fp_points = sum([s['fp_count'] for s in slide_details])
    avg_fp_points = total_fp_points / len(slide_details) if slide_details else 0
    print(f"\nFalse positive point statistics:")
    print(f"Total FP points: {total_fp_points}")
    print(f"Average FP points per slide: {avg_fp_points:.2f}")


# if __name__ == "__main__":
#     # Usage example
#     overall_power, fp_rate, slide_details, file_count = calculate_power(
#         folder_path=r"C:\Users\pmgl903\Desktop\kgg\kggsum3\simulation\result_scdrs_32456",
#         disease_regions=[(20, 29, 30, 39), (40,49,30,39), (60,69,50,59)], # Two disease regions
#         inner_regions=[(14, 35, 24, 45), (4,55,24,45), (54,75,44,65)] ,     # Two inner regions: X-axis 0-99, Y-axis 14-35 and 54-75
#         power_threshold=0.8,  # Slide power is 1 when 80% points are significant
#         sig_threshold=0.05
#     )
    
#     # Print detailed results
#     print_detailed_results(slide_details, overall_power, fp_rate, file_count)