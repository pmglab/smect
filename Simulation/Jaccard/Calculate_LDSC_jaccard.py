import pandas as pd
import numpy as np
import os
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def load_and_process_ldsc_data(correct_dir, wrong_dir, alpha=0.05):
    """
    Load and process data with correct labels and wrong labels
    
    Parameters:
    -----------
    correct_dir : str
        Path to the folder containing correct label results
    wrong_dir : str
        Path to the folder containing wrong label results
    alpha : float
        Significance level threshold
    
    Returns:
    --------
    dict : Dictionary containing processed results
    """
    
    # Get file lists
    correct_files = [f for f in os.listdir(correct_dir) if f.endswith('.gz')]
    wrong_files = [f for f in os.listdir(wrong_dir) if f.endswith('.gz')]
    
    # Ensure file correspondence
    correct_files.sort()
    wrong_files.sort()
    
    results = {
        'correct_significant_spots': [],
        'wrong_significant_spots': [],
        'jaccard_similarities': [],
        'overlap_ratios': [],
        'file_pairs': []
    }
    
    for correct_file, wrong_file in zip(correct_files, wrong_files):
        # Read correct label data
        correct_path = os.path.join(correct_dir, correct_file)
        correct_df = pd.read_csv(correct_path)
        
        # Read wrong label data
        wrong_path = os.path.join(wrong_dir, wrong_file)
        wrong_df = pd.read_csv(wrong_path)
        
        # Multiple testing correction (FDR correction)
        _, correct_pvals_corrected, _, _ = multipletests(
            correct_df['p'], alpha=alpha, method='fdr_bh'
        )
        
        _, wrong_pvals_corrected, _, _ = multipletests(
            wrong_df['p'], alpha=alpha, method='fdr_bh'
        )
        
        # Add corrected p-values
        correct_df['p_corrected'] = correct_pvals_corrected
        wrong_df['p_corrected'] = wrong_pvals_corrected
        
        # Identify significant spots
        correct_significant = set(correct_df[correct_df['p_corrected'] < alpha]['spot'])
        wrong_significant = set(wrong_df[wrong_df['p_corrected'] < alpha]['spot'])
        
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
        results['correct_significant_spots'].append(correct_significant)
        results['wrong_significant_spots'].append(wrong_significant)
        results['jaccard_similarities'].append(jaccard)
        results['overlap_ratios'].append(overlap_ratio)
        results['file_pairs'].append((correct_file, wrong_file))
    
    return results

def calculate_comprehensive_metrics(results):
    """
    Calculate comprehensive evaluation metrics
    """
    metrics = {}
    
    # Basic statistics
    metrics['mean_jaccard'] = np.mean(results['jaccard_similarities'])
    metrics['std_jaccard'] = np.std(results['jaccard_similarities'])
    metrics['mean_overlap_ratio'] = np.mean(results['overlap_ratios'])
    metrics['std_overlap_ratio'] = np.std(results['overlap_ratios'])
    
    # Stability metrics
    metrics['jaccard_cv'] = metrics['std_jaccard'] / metrics['mean_jaccard'] if metrics['mean_jaccard'] > 0 else np.inf
    metrics['overlap_ratio_cv'] = metrics['std_overlap_ratio'] / metrics['mean_overlap_ratio'] if metrics['mean_overlap_ratio'] > 0 else np.inf
    
    return metrics

def plot_comparison_results(results, metrics):
    """
    Plot visualization charts for comparison results
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Jaccard similarity coefficient distribution
    axes[0, 0].hist(results['jaccard_similarities'], bins=20, alpha=0.7, color='skyblue')
    axes[0, 0].axvline(metrics['mean_jaccard'], color='red', linestyle='--', 
                      label=f'Mean jaccard: {metrics["mean_jaccard"]:.3f}')
    axes[0, 0].set_xlabel('Jaccard similarity coefficient')
    axes[0, 0].set_ylabel('frequency')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Overlap ratio distribution
    axes[0, 1].hist(results['overlap_ratios'], bins=20, alpha=0.7, color='lightgreen')
    axes[0, 1].axvline(metrics['mean_overlap_ratio'], color='red', linestyle='--',
                      label=f'Mean overlap ratio: {metrics["mean_overlap_ratio"]:.3f}')
    axes[0, 1].set_xlabel('overlap ratios')
    axes[0, 1].set_ylabel('frequency')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Relationship between two metrics
    axes[1, 0].scatter(results['jaccard_similarities'], results['overlap_ratios'], alpha=0.6)
    axes[1, 0].set_xlabel('Jaccard similarity coefficient')
    axes[1, 0].set_ylabel('overlap ratios')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Metric variation across file pairs
    file_indices = range(len(results['jaccard_similarities']))
    axes[1, 1].plot(file_indices, results['jaccard_similarities'], 'o-', label='Jaccard similarity coefficient')
    axes[1, 1].plot(file_indices, results['overlap_ratios'], 's-', label='overlap ratios')
    axes[1, 1].set_xlabel('File pair index')
    axes[1, 1].set_ylabel('Metric value')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('C:/Users/pmgl903/Desktop/研究生/论文/Fig/simulation/fig3.png', dpi=800)
    plt.show()

def print_detailed_report(results, metrics):
    """
    Print detailed evaluation report
    """
    print("="*60)
    print("             Anti-interference Capability Evaluation Report")
    print("="*60)
    
    print(f"\n1. Basic Statistics:")
    print(f"   Number of file pairs: {len(results['file_pairs'])}")
    print(f"   Jaccard similarity coefficient - Mean: {metrics['mean_jaccard']:.4f} ± {metrics['std_jaccard']:.4f}")
    print(f"   Overlap ratio - Mean: {metrics['mean_overlap_ratio']:.4f} ± {metrics['std_overlap_ratio']:.4f}")
    
    print(f"\n2. Stability Metrics:")
    print(f"   Jaccard coefficient of variation: {metrics['jaccard_cv']:.4f}")
    print(f"   Overlap ratio coefficient of variation: {metrics['overlap_ratio_cv']:.4f}")
    
    print(f"\n3. Detailed File Pair Results:")
    for i, (jaccard, overlap, (correct_file, wrong_file)) in enumerate(
        zip(results['jaccard_similarities'], results['overlap_ratios'], results['file_pairs'])
    ):
        correct_count = len(results['correct_significant_spots'][i])
        wrong_count = len(results['wrong_significant_spots'][i])
        overlap_count = len(results['correct_significant_spots'][i] & results['wrong_significant_spots'][i])
        
        print(f"   File pair {i:2d}: Jaccard={jaccard:.4f}, Overlap ratio={overlap:.4f}")
        print(f"           Correct significant: {correct_count:3d}, Wrong significant: {wrong_count:3d}, Overlap: {overlap_count:3d}")