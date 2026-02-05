import pandas as pd
import numpy as np
import os
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

def load_enrichment_dese_data(correct_dir, wrong_dir, alpha=0.05):
    """
    Load enrichment analysis result data
    
    Parameters:
    -----------
    correct_dir : str
        Path to the correct label results folder
    wrong_dir : str
        Path to the wrong label results folder
    alpha : float
        Significance threshold
    
    Returns:
    --------
    dict : Dictionary containing processed results
    """
    
    # Get file list
    correct_files = sorted(glob.glob(os.path.join(correct_dir, "*.enrichment.xls")))
    wrong_files = sorted(glob.glob(os.path.join(wrong_dir, "*.enrichment.xls")))
    
    print(f"Found correct label files: {len(correct_files)}")
    print(f"Found wrong label files: {len(wrong_files)}")
    
    # Extract file numbers for matching
    def extract_file_number(filepath):
        filename = os.path.basename(filepath)
        # Extract numbers from filename, e.g., extract 0 from "gene.spatial.expression.txt.0.enrichment.xls"
        parts = filename.split('.')
        for part in parts:
            if part.isdigit():
                return int(part)
        return -1
    
    # Organize files by file number
    correct_files_dict = {extract_file_number(f): f for f in correct_files}
    wrong_files_dict = {extract_file_number(f): f for f in wrong_files}
    
    # Find common file numbers
    common_numbers = set(correct_files_dict.keys()) & set(wrong_files_dict.keys())
    common_numbers = sorted([n for n in common_numbers if n >= 0])
    
    print(f"Common file numbers: {common_numbers}")
    
    results = {
        'correct_significant_conditions': [],
        'wrong_significant_conditions': [],
        'jaccard_similarities': [],
        'overlap_ratios': [],
        'file_pairs': [],
        'all_conditions': set(),
        'correct_p_values': [],
        'wrong_p_values': []
    }
    
    for file_num in common_numbers:
        correct_file = correct_files_dict[file_num]
        wrong_file = wrong_files_dict[file_num]
        
        try:
            # Read correct label data
            correct_df = pd.read_excel(correct_file)
            # Read wrong label data
            wrong_df = pd.read_excel(wrong_file)
            
            print(f"\nProcessing file number {file_num}:")
            print(f"  Correct file: {os.path.basename(correct_file)}, row count: {len(correct_df)}")
            print(f"  Wrong file: {os.path.basename(wrong_file)}, row count: {len(wrong_df)}")
            
            # Check column names
            correct_columns = correct_df.columns.tolist()
            wrong_columns = wrong_df.columns.tolist()
            
            print(f"  Correct file column names: {correct_columns}")
            print(f"  Wrong file column names: {wrong_columns}")
            
            # Based on image information, column names should be "Condition" and "Adjusted(p)"
            condition_col = "Condition"
            pval_col = "Adjusted(p)"
            
            # Verify if column names exist
            if condition_col not in correct_columns or pval_col not in correct_columns:
                print(f"  Warning: Column names don't match, trying to use first and fourth columns")
                if len(correct_columns) >= 4:
                    condition_col = correct_columns[0]
                    pval_col = correct_columns[3]
                else:
                    print(f"  Error: File has insufficient columns, skipping")
                    continue
            
            print(f"  Using condition column: '{condition_col}', p-value column: '{pval_col}'")
            
            # Extract significant conditions (using threshold 0.05)
            correct_significant = set(correct_df[correct_df[pval_col] < alpha][condition_col].astype(str).unique())
            wrong_significant = set(wrong_df[wrong_df[pval_col] < alpha][condition_col].astype(str).unique())
            
            # Record p-value distribution
            results['correct_p_values'].extend(correct_df[pval_col].tolist())
            results['wrong_p_values'].extend(wrong_df[pval_col].tolist())
            
            # Record all conditions
            results['all_conditions'].update(correct_df[condition_col].astype(str).unique())
            results['all_conditions'].update(wrong_df[condition_col].astype(str).unique())
            
            # Calculate overlap metrics
            intersection = correct_significant & wrong_significant
            union = correct_significant | wrong_significant
            
            # Jaccard similarity coefficient
            if len(union) > 0:
                jaccard = len(intersection) / len(union)
            else:
                jaccard = 0
            
            # Overlap ratio (based on correct labels)
            if len(correct_significant) > 0:
                overlap_ratio = len(intersection) / len(correct_significant)
            else:
                overlap_ratio = 0
            
            # Store results
            results['correct_significant_conditions'].append(correct_significant)
            results['wrong_significant_conditions'].append(wrong_significant)
            results['jaccard_similarities'].append(jaccard)
            results['overlap_ratios'].append(overlap_ratio)
            results['file_pairs'].append((correct_file, wrong_file))
            
            print(f"  Correct significant conditions: {len(correct_significant)}")
            print(f"  Wrong significant conditions: {len(wrong_significant)}")
            print(f"  Overlapping conditions: {len(intersection)}")
            print(f"  Jaccard similarity coefficient: {jaccard:.4f}")
            print(f"  Overlap ratio: {overlap_ratio:.4f}")
            
        except Exception as e:
            print(f"  Error processing file {file_num}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    return results

def calculate_comprehensive_metrics(results):
    """
    Calculate comprehensive evaluation metrics
    """
    metrics = {}
    
    if len(results['jaccard_similarities']) == 0:
        print("Warning: No file pairs were successfully processed")
        return metrics
    
    # Basic statistics
    metrics['mean_jaccard'] = np.mean(results['jaccard_similarities'])
    metrics['std_jaccard'] = np.std(results['jaccard_similarities'])
    metrics['mean_overlap_ratio'] = np.mean(results['overlap_ratios'])
    metrics['std_overlap_ratio'] = np.std(results['overlap_ratios'])
    
    # Stability metrics
    metrics['jaccard_cv'] = metrics['std_jaccard'] / metrics['mean_jaccard'] if metrics['mean_jaccard'] > 0 else np.inf
    metrics['overlap_ratio_cv'] = metrics['std_overlap_ratio'] / metrics['mean_overlap_ratio'] if metrics['mean_overlap_ratio'] > 0 else np.inf
    
    # File pair count
    metrics['file_pairs_count'] = len(results['file_pairs'])
    
    # p-value statistics
    if results['correct_p_values']:
        metrics['correct_p_mean'] = np.mean(results['correct_p_values'])
        metrics['correct_p_std'] = np.std(results['correct_p_values'])
        metrics['correct_p_significant'] = sum(1 for p in results['correct_p_values'] if p < 0.05) / len(results['correct_p_values'])
    
    if results['wrong_p_values']:
        metrics['wrong_p_mean'] = np.mean(results['wrong_p_values'])
        metrics['wrong_p_std'] = np.std(results['wrong_p_values'])
        metrics['wrong_p_significant'] = sum(1 for p in results['wrong_p_values'] if p < 0.05) / len(results['wrong_p_values'])
    
    return metrics

def plot_comparison_results(results, metrics):
    """
    Plot visualization charts for comparison results
    """
    if len(results['jaccard_similarities']) == 0:
        print("No data available for plotting")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Jaccard similarity coefficient distribution
    axes[0, 0].hist(results['jaccard_similarities'], bins=20, alpha=0.7, color='skyblue')
    if 'mean_jaccard' in metrics:
        axes[0, 0].axvline(metrics['mean_jaccard'], color='red', linestyle='--', 
                          label=f'Mean jaccard: {metrics["mean_jaccard"]:.3f}')
    axes[0, 0].set_xlabel('Jaccard similarity coefficient')
    axes[0, 0].set_ylabel('frequency')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('C:/Users/pmgl903/Desktop/研究生/论文/Fig/simulation/fig3.png',dpi=800)
    plt.show()

def print_detailed_report(results, metrics):
    """
    Print detailed evaluation report
    """
    print("="*60)
    print("             Anti-interference Capability Evaluation Report")
    print("="*60)
    
    if len(results['file_pairs']) == 0:
        print("No available file pairs for analysis")
        return
    
    print(f"\n1. Basic Statistics:")
    print(f"   Successfully processed file pairs: {len(results['file_pairs'])}")
    if metrics:
        print(f"   Jaccard similarity coefficient - Mean: {metrics['mean_jaccard']:.4f} ± {metrics['std_jaccard']:.4f}")
        print(f"   Overlap ratio - Mean: {metrics['mean_overlap_ratio']:.4f} ± {metrics['std_overlap_ratio']:.4f}")
    
    print(f"\n2. Stability Metrics:")
    if metrics:
        print(f"   Jaccard coefficient of variation: {metrics['jaccard_cv']:.4f}")
        print(f"   Overlap ratio coefficient of variation: {metrics['overlap_ratio_cv']:.4f}")
    
    print(f"\n3. P-value Statistics:")
    if metrics:
        if 'correct_p_mean' in metrics:
            print(f"   Correct label p-values - Mean: {metrics['correct_p_mean']:.4f} ± {metrics['correct_p_std']:.4f}")
            print(f"   Correct label significant proportion: {metrics['correct_p_significant']:.4f}")
        if 'wrong_p_mean' in metrics:
            print(f"   Wrong label p-values - Mean: {metrics['wrong_p_mean']:.4f} ± {metrics['wrong_p_std']:.4f}")
            print(f"   Wrong label significant proportion: {metrics['wrong_p_significant']:.4f}")
    
    print(f"\n4. Detailed File Pair Results:")
    for i, (jaccard, overlap, (correct_file, wrong_file)) in enumerate(
        zip(results['jaccard_similarities'], results['overlap_ratios'], results['file_pairs'])
    ):
        correct_count = len(results['correct_significant_conditions'][i])
        wrong_count = len(results['wrong_significant_conditions'][i])
        overlap_count = len(results['correct_significant_conditions'][i] & results['wrong_significant_conditions'][i])
        
        print(f"   File pair {i:2d}: Jaccard={jaccard:.4f}, Overlap ratio={overlap:.4f}")
        print(f"           Correct significant: {correct_count:3d}, Wrong significant: {wrong_count:3d}, Overlap: {overlap_count:3d}")

def analyze_condition_patterns(results):
    """
    Analyze condition patterns to find common and unique conditions
    """
    if len(results['correct_significant_conditions']) == 0:
        return
    
    print(f"\n5. Condition Pattern Analysis:")
    
    # Find conditions that are common across all file pairs
    if len(results['correct_significant_conditions']) > 0:
        common_conditions = set.intersection(*[correct & wrong for correct, wrong in 
                                             zip(results['correct_significant_conditions'], 
                                                 results['wrong_significant_conditions'])])
        
        print(f"   Conditions significant in all file pairs: {len(common_conditions)}")
        if common_conditions:
            print(f"   Examples: {sorted(common_conditions)[:10]}")  # Show only first 10
    
    # Analyze condition stability
    condition_stability = defaultdict(int)
    for correct_conditions in results['correct_significant_conditions']:
        for condition in correct_conditions:
            condition_stability[condition] += 1
    
    if len(results['correct_significant_conditions']) > 0:
        stable_conditions = [cond for cond, count in condition_stability.items() 
                            if count >= len(results['correct_significant_conditions']) * 0.8]
        
        print(f"   Conditions significant in over 80% of files: {len(stable_conditions)}")