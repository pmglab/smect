import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import gzip
import pathlib
from typing import List, Union, Optional
import scanpy as sc

class calculate_OR:
    """
    A class for calculating Odds Ratios (OR), 95% Confidence Intervals, and p-values
    Supports analysis of DESE, LDSC, and scDRS data formats
    """
    
    def __init__(self, adata):
        """
        Initialize the calculate_OR class
        
        Parameters:
            adata: AnnData object containing single-cell data
        """
        self.adata = adata
        self.results = {}
        print(f"calculate_OR class initialized successfully, containing {adata.n_obs} samples")
    
    def process_enrichment_DESE_files(self, folder_path, col_indices=[0, 3]):
        """
        Batch process DESE enrichment analysis files and merge into obs
        
        Parameters:
            folder_path: Path to folder containing .xls files
            col_indices: Column indices to retain, default [0,3] (1st and 4th columns)
        """
        # Get all .xls file paths
        file_list = [os.path.join(folder_path, f) 
                    for f in os.listdir(folder_path) 
                    if f.endswith(".xls")]
        
        processed_count = 0
        for file_path in file_list:
            try:
                # Read data
                df = pd.read_excel(file_path, engine="xlrd")
                
                # Extract filename features as suffix
                filename = Path(file_path).name
                suffix = filename.split("_gss_")[-1].split(".")[0]
                
                # Process dataframe
                new_df = df.iloc[:, col_indices]
                new_df = new_df.set_index(new_df.columns[0])
                new_df = new_df.drop(new_df.index[-2:])
                new_df.index.name = 'cell_name'
                
                # Rename column with suffix
                original_col = new_df.columns[0]
                new_colname = f"{original_col}_{suffix}"
                new_df = new_df.rename(columns={original_col: new_colname})
                
                # Merge into obs
                valid_cells = new_df.index.intersection(self.adata.obs.index)
                self.adata.obs = self.adata.obs.join(
                    new_df.loc[valid_cells], 
                    how='left'
                )
                processed_count += 1
                
            except Exception as e:
                print(f"Error processing file {file_path}: {str(e)}")
                continue
        
        print(f"Successfully merged {processed_count} DESE files into obs matrix")
        return self.adata
    
    def process_enrichment_LDSC_files(self, folder_path, col_indices=[0, 4]):
        """
        Batch process LDSC enrichment analysis files and merge into obs
        
        Parameters:
            folder_path: Path to folder containing .csv files
            col_indices: Column indices to retain, default [0,4] (1st and 5th columns)
        """
        # Get all .csv file paths
        file_list = [os.path.join(folder_path, f) 
                    for f in os.listdir(folder_path) 
                    if f.endswith(".csv") and "gsmap" in f]
        
        processed_count = 0
        for file_path in file_list:
            try:
                # Read CSV file
                df = pd.read_csv(file_path)
                
                # Extract filename features as suffix
                filename = Path(file_path).name
                suffix = filename.split("_gsmap_")[-1].split(".")[0]
                
                # Process dataframe
                new_df = df.iloc[:, col_indices]
                new_df = new_df.set_index(new_df.columns[0])
                new_df.index.name = 'cell_name'
                
                # Rename column with disease identifier suffix
                original_col = new_df.columns[0]
                new_colname = f"{original_col}_{suffix}"
                new_df = new_df.rename(columns={original_col: new_colname})
                
                # Merge into obs
                valid_cells = new_df.index.intersection(self.adata.obs.index)
                self.adata.obs = self.adata.obs.join(
                    new_df.loc[valid_cells], 
                    how='left'
                )
                processed_count += 1
                
            except Exception as e:
                print(f"Error processing file {file_path}: {str(e)}")
                continue
        
        print(f"Successfully merged {processed_count} LDSC files into obs matrix")
        return self.adata
    
    def merge_scDRS_scores(self, input_dir, files_to_process=None, save_path=None, p_value_col='mc_pval'):
        """
        Merge scDRS score files into the obs of AnnData object
        
        Parameters:
            input_dir: Directory path containing scDRS score files
            files_to_process: List of score filenames to process. If None, all .score.gz files will be processed
            save_path: Optional path to save the merged AnnData object
            p_value_col: Column name to use for p-values. Options: 'mc_pval', 'pval', 'nlog10_pval'
            
        Returns:
            Updated AnnData object with scDRS p-values added to obs
            
        Example:
            >>> or_calculator = calculate_OR(mmdata)
            >>> updated_adata = or_calculator.merge_scDRS_scores("path/to/scores")
        """
        # Check if input directory exists
        input_dir = Path(input_dir)
        if not input_dir.exists():
            print(f"Error: Directory '{input_dir}' does not exist!")
            print(f"Please check the path and make sure it's correct.")
            return self.adata
        
        if not input_dir.is_dir():
            print(f"Error: '{input_dir}' is not a directory!")
            return self.adata
        
        print(f"Processing scDRS files from directory: {input_dir}")
        print(f"Directory exists: {input_dir.exists()}")
        
        # If files_to_process is not provided, get all .score.gz files in the directory
        if files_to_process is None:
            try:
                all_files = os.listdir(input_dir)
                files_to_process = [f for f in all_files 
                                  if f.endswith('.score.gz') and os.path.isfile(input_dir / f)]
                
                if not files_to_process:
                    print(f"Warning: No .score.gz files found in '{input_dir}'!")
                    print(f"Available files: {all_files[:10]}")  # Show first 10 files for debugging
                    return self.adata
                
                print(f"Found {len(files_to_process)} .score.gz files to process")
            except Exception as e:
                print(f"Error listing files in directory: {str(e)}")
                return self.adata
        elif isinstance(files_to_process, list) and len(files_to_process) == 0:
            print("Warning: files_to_process is an empty list, no files to process")
            return self.adata
        
        # Validate p_value_col parameter based on image data
        valid_p_value_cols = ['mc_pval', 'pval', 'nlog10_pval']
        if p_value_col not in valid_p_value_cols:
            print(f"Warning: p_value_col must be one of {valid_p_value_cols}, using 'mc_pval' as default")
            p_value_col = 'mc_pval'
        
        print(f"Using p-value column: {p_value_col}")
        print(f"Processing {len(files_to_process)} files")
        
        # Store the list of successfully processed files
        processed_files = []
        failed_files = []
        
        for file_idx, file_name in enumerate(files_to_process, 1):
            file_path = input_dir / file_name
            
            if not file_path.exists():
                print(f"\n[{file_idx}/{len(files_to_process)}] Error: File '{file_name}' does not exist!")
                failed_files.append(file_name)
                continue
                
            # Extract trait name from filename (remove .score.gz)
            trait = file_name.replace('.score.gz', '')
            col_name = f"p_{trait}"
            
            print(f"\n[{file_idx}/{len(files_to_process)}] Processing: {file_name}")
            print(f"  Trait: {trait}")
            print(f"  Column name: {col_name}")
            print(f"  File path: {file_path}")
            print(f"  File exists: {file_path.exists()}")
            
            try:
                # Read the gzipped score file
                print(f"  Reading file...")
                with gzip.open(file_path, 'rt') as f:
                    # First, read the first few lines to understand the structure
                    preview_lines = []
                    for _ in range(5):
                        line = f.readline()
                        if not line:
                            break
                        preview_lines.append(line.strip())
                    
                    # Reset file pointer
                    f.seek(0)
                    
                    # Debug: show preview of file
                    print(f"  First few lines of file:")
                    for i, line in enumerate(preview_lines[:3]):
                        print(f"    Line {i+1}: {line[:100]}...")
                    
                    # Check if it's tab-separated
                    if '\t' in preview_lines[0]:
                        sep = '\t'
                        print(f"  Detected tab-separated format")
                    else:
                        # Try to detect separator
                        if ',' in preview_lines[0]:
                            sep = ','
                            print(f"  Detected comma-separated format")
                        else:
                            # Use whitespace as separator
                            sep = r'\s+'
                            print(f"  Detected whitespace-separated format")
                    
                    # Read the data
                    df = pd.read_csv(f, sep=sep, engine='python')
                    
                    print(f"  File shape: {df.shape}")
                    print(f"  Columns: {list(df.columns)}")
                    
                    # Check if we have the expected columns
                    if p_value_col not in df.columns:
                        print(f"  Warning: Column '{p_value_col}' not found in file!")
                        print(f"  Available columns: {list(df.columns)}")
                        
                        # Try to find alternative p-value columns
                        alternative_cols = [c for c in valid_p_value_cols if c in df.columns]
                        if alternative_cols:
                            p_value_col_to_use = alternative_cols[0]
                            print(f"  Using alternative column: '{p_value_col_to_use}'")
                        else:
                            print(f"  Skipping file - no suitable p-value columns found")
                            failed_files.append(file_name)
                            continue
                    else:
                        p_value_col_to_use = p_value_col
                    
                    # Handle the index column
                    # Based on the image, the first column doesn't have a header name
                    # It contains values like '190_194', '190_195', etc.
                    if df.columns[0] == 'Unnamed: 0':
                        # Set the first column as index
                        df.set_index('Unnamed: 0', inplace=True)
                        df.index.name = 'cell_name'
                        print(f"  Set 'Unnamed: 0' column as index")
                    elif df.index.name is None and len(df.columns) > 0:
                        # Check if the first column looks like a cell/sample identifier
                        first_col = df.columns[0]
                        if any(char in str(first_col) for char in ['_', '-']) or df[first_col].astype(str).str.contains('_').any():
                            df.set_index(first_col, inplace=True)
                            df.index.name = 'cell_name'
                            print(f"  Set '{first_col}' column as index")
                    
                    # If still no index name, set a default
                    if df.index.name is None:
                        df.index.name = 'cell_name'
                    
                    print(f"  Index name: {df.index.name}")
                    print(f"  Sample index examples: {df.index[:5].tolist()}")
                    
                    # Convert p-values if using nlog10_pval (negative log10 p-values)
                    if p_value_col_to_use == 'nlog10_pval':
                        print(f"  Converting nlog10_pval to linear p-values")
                        # nlog10_pval = -log10(p), so p = 10^(-nlog10_pval)
                        p_values = 10 ** (-df[p_value_col_to_use])
                    else:
                        p_values = df[p_value_col_to_use]
                    
                    # Check sample ID matching
                    common_samples = df.index.intersection(self.adata.obs.index)
                    print(f"  Matching samples: {len(common_samples)} out of {len(df)}")
                    
                    if len(common_samples) == 0:
                        print(f"  Warning: No matching samples found!")
                        print(f"  First 5 file indices: {df.index[:5].tolist()}")
                        print(f"  First 5 adata indices: {self.adata.obs.index[:5].tolist()}")
                        
                        # Try to see if there's a pattern mismatch
                        # Sometimes indices might have different prefixes/suffixes
                        adata_samples_str = ' '.join(self.adata.obs.index[:5].astype(str))
                        file_samples_str = ' '.join(df.index[:5].astype(str))
                        print(f"  String representation comparison:")
                        print(f"    Adata: {adata_samples_str}")
                        print(f"    File:  {file_samples_str}")
                        
                        failed_files.append(file_name)
                        continue
                    
                    # Check for duplicate indices in the file
                    if df.index.duplicated().any():
                        print(f"  Warning: Duplicate indices found in file, keeping first occurrence")
                        df = df[~df.index.duplicated(keep='first')]
                        # Update common_samples after deduplication
                        common_samples = df.index.intersection(self.adata.obs.index)
                    
                    # Create or update the column
                    if col_name in self.adata.obs.columns:
                        print(f"  Column {col_name} already exists, updating values")
                        # Update only the matching samples
                        self.adata.obs.loc[common_samples, col_name] = p_values.loc[common_samples].values
                    else:
                        print(f"  Creating new column: {col_name}")
                        # Initialize with NaN
                        self.adata.obs[col_name] = np.nan
                        # Fill in values for matching samples
                        self.adata.obs.loc[common_samples, col_name] = p_values.loc[common_samples].values
                    
                    # Report statistics
                    non_null_count = self.adata.obs[col_name].notna().sum()
                    if non_null_count > 0:
                        col_mean = self.adata.obs[col_name].mean()
                        col_min = self.adata.obs[col_name].min()
                        col_max = self.adata.obs[col_name].max()
                        print(f"  Added {col_name}:")
                        print(f"    Non-NaN values: {non_null_count}")
                        print(f"    Mean p-value: {col_mean:.6f}")
                        print(f"    Min p-value: {col_min:.6f}")
                        print(f"    Max p-value: {col_max:.6f}")
                    else:
                        print(f"  Warning: No values were added to {col_name}!")
                    
                    processed_files.append(file_name)
                    
            except Exception as e:
                print(f"  Error processing file {file_name}: {str(e)}")
                import traceback
                traceback.print_exc()
                failed_files.append(file_name)
                continue
        
        # Report final status
        print(f"\n{'='*80}")
        print("SCDRS SCORES MERGING SUMMARY")
        print(f"{'='*80}")
        print(f"Total files processed: {len(processed_files)}")
        print(f"Successfully processed: {len(processed_files)}")
        print(f"Failed to process: {len(failed_files)}")
        
        if processed_files:
            print(f"\nSuccessfully processed files:")
            for f in processed_files:
                print(f"  ✓ {f}")
        
        if failed_files:
            print(f"\nFailed files:")
            for f in failed_files:
                print(f"  ✗ {f}")
        
        # Count p-value columns
        p_cols = [col for col in self.adata.obs.columns if col.startswith('p_')]
        print(f"\nTotal p_* columns in obs: {len(p_cols)}")
        if p_cols:
            print(f"p-value columns: {p_cols}")
            
            # Show summary statistics for p-value columns
            print(f"\nSummary statistics for p-value columns:")
            for col in sorted(p_cols):
                non_null = self.adata.obs[col].notna().sum()
                if non_null > 0:
                    mean_val = self.adata.obs[col].mean()
                    min_val = self.adata.obs[col].min()
                    max_val = self.adata.obs[col].max()
                    print(f"  {col}: {non_null} values, mean={mean_val:.6f}, min={min_val:.6f}, max={max_val:.6f}")
                else:
                    print(f"  {col}: No values (all NaN)")
        
        # Save updated AnnData object if requested
        if save_path:
            try:
                save_path = Path(save_path)
                self.adata.write(save_path)
                print(f"\nUpdated AnnData object saved to: {save_path}")
            except Exception as e:
                print(f"\nError saving AnnData object: {str(e)}")
        
        return self.adata
    
    def batch_multiple_test_correction(self, p_prefix='p_', method='fdr_bh'):
        """
        Perform batch multiple testing correction on p-value columns with specified prefix
        
        Parameters:
            p_prefix: Prefix of p-value columns
            method: Multiple testing correction method, default 'fdr_bh'
        """
        from statsmodels.stats.multitest import multipletests
        import numpy as np
        
        # Find all column names starting with the specified prefix
        p_cols = [col for col in self.adata.obs.columns if col.startswith(p_prefix)]
        
        if not p_cols:
            print(f"No columns found starting with '{p_prefix}'")
            return self.adata
        
        print(f"Found {len(p_cols)} p-value columns requiring correction")
        
        for col in p_cols:
            original_p_values = self.adata.obs[col].values
            
            # === Key modification: Create mask for non-NaN values ===
            non_nan_mask = ~np.isnan(original_p_values)
            valid_p_values = original_p_values[non_nan_mask]  # Extract only valid p-values
            
            # Check if there are valid p-values to process
            if len(valid_p_values) == 0:
                print(f"Column '{col}' contains no valid p-values (all NaN). Skipping.")
                # Option to create a column of all NaN or skip directly
                new_col = f'p_adj_{col.split("_")[-1]}'
                self.adata.obs[new_col] = np.nan  # Explicitly create a column of NaN
                continue
            
            # Perform multiple testing correction only on valid p-values
            reject, p_adj = multipletests(valid_p_values, method=method)[:2]
            
            # Create new column name
            disease_abbr = col.split('_')[-1]
            new_col = f'p_adj_{disease_abbr}'
            
            # === Key modification: Create all-NaN array and fill back valid positions ===
            adjusted_values = np.full_like(original_p_values, np.nan, dtype=float)
            adjusted_values[non_nan_mask] = p_adj
            
            # Add the corrected p-value array (including NaN) to the DataFrame
            self.adata.obs[new_col] = adjusted_values
            
            # Statistical information (calculated based on valid values)
            print(f"Processed column {col} -> {new_col}")
            print(f"  Valid p-values: {len(valid_p_values)}")
            print(f"  Significant results (after correction): {sum(reject)}")
        
        return self.adata
    
    def _calculate_merged_or(self, target_regions, trait_col='trait_related'):
        """
        Calculate OR values, 95% CI, and p-values for merged brain regions (internal method)
        
        Parameters:
            target_regions: List of target brain regions
            trait_col: Trait column name
        """
        # Generate merged region mask
        merged_mask = self.adata.obs['annotation'].isin(target_regions)
        
        # Build contingency table (add pseudo-counts for zero values)
        raw_table = pd.crosstab(
            merged_mask,
            self.adata.obs[trait_col].astype(bool)
        ).reindex(index=[True, False], columns=[True, False]).fillna(0)
        
        # Add pseudo-counts
        cont_table = raw_table + 0.5
        
        # Calculate OR and CI
        a, b = cont_table.loc[True]
        c, d = cont_table.loc[False]
        
        try:
            OR = (a * d) / (b * c)
            se = np.sqrt(1/a + 1/b + 1/c + 1/d)
            ci_lower = np.exp(np.log(OR) - 1.96*se)
            ci_upper = np.exp(np.log(OR) + 1.96*se)
        except ZeroDivisionError:
            OR = np.nan
            ci_lower, ci_upper = np.nan, np.nan
        
        # Chi-square test
        _, p, _, _ = chi2_contingency(raw_table.values, correction=False)
        
        return OR, (ci_lower, ci_upper), p
    
    def batch_analyze_traits(self, analysis_type="DESE", target_regions_list=None, 
                           output_csv="OR_results.csv", p_value_threshold=0.05):
        """
        Main analysis function supporting three analysis types
        
        Parameters:
            analysis_type: Analysis type, options: "DESE", "LDSC", "scDRS"
            target_regions_list: List of brain regions to analyze
            output_csv: Path to save results
            p_value_threshold: p-value threshold, default 0.05
            
        Returns:
            Two DataFrames: formatted_df (matching image format) and detailed_df (detailed results)
        """
        # Set default brain region list
        if target_regions_list is None:
            unique_regions = self.adata.obs['annotation'].dropna().unique().tolist()
            target_regions_list = [[region] for region in unique_regions]
        
        # Apply multiple testing correction for LDSC and scDRS, but not for DESE
        if analysis_type == "LDSC":
            print("Applying multiple testing correction for LDSC analysis...")
            self.batch_multiple_test_correction(p_prefix='p_gsmap_E1S1_10_', method='fdr_bh')
        elif analysis_type == "scDRS":
            print("Applying multiple testing correction for scDRS analysis...")
            self.batch_multiple_test_correction(p_prefix='p_', method='fdr_bh')
        elif analysis_type == "DESE":
            print("DESE analysis: No multiple testing correction applied.")
        else:
            raise ValueError("analysis_type must be 'DESE', 'LDSC', or 'scDRS'")
        
        # Determine trait column prefix based on analysis type
        if analysis_type == "DESE":
            p_prefix = "Adjusted(p)_"
            trait_prefix = "trait_"
        elif analysis_type == "LDSC":
            p_prefix = "p_adj_"
            trait_prefix = "trait_"
        elif analysis_type == "scDRS":
            p_prefix = "p_adj_"
            trait_prefix = "trait_"
        
        # Automatically get list of traits
        traits = [col.split("_")[-1] 
                 for col in self.adata.obs.columns 
                 if col.startswith(p_prefix)]
        
        if not traits:
            print(f"No trait columns found starting with '{p_prefix}'")
            return None, None
        
        results = []
        
        # Iterate through each trait
        for trait in traits:
            p_col = f'{p_prefix}{trait}'
            trait_col = f'{trait_prefix}{trait}'
            
            # Create temporary trait column
            self.adata.obs[trait_col] = (self.adata.obs[p_col] < p_value_threshold)
            
            # Iterate through each brain region combination
            for regions in target_regions_list:
                region_name = "+".join(regions)
                
                # Calculate OR
                try:
                    OR, (ci_low, ci_high), p = self._calculate_merged_or(regions, trait_col)
                except Exception as e:
                    print(f"Error in {trait}-{region_name}: {str(e)}")
                    OR, ci_low, ci_high, p = np.nan, np.nan, np.nan, np.nan
                    
                # Store results
                results.append({
                    'Region': region_name,
                    'Trait': trait,
                    'OR': OR,
                    'CI_lower': ci_low,
                    'CI_upper': ci_high,
                    'p_value': p,
                    'Analysis_Type': analysis_type
                })
            
            # Clean up temporary column
            del self.adata.obs[trait_col]
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Create two output formats:
        # 1. Detailed format (long format)
        detailed_df = results_df.copy()
        
        # 2. Format matching the image (wide format with OR, CI_lower, CI_upper columns for each trait)
        formatted_results = []
        
        # Get unique regions
        unique_regions = sorted(results_df['Region'].unique())
        
        for region in unique_regions:
            region_data = {'Region': region}
            
            # Get all traits for this region
            region_traits = results_df[results_df['Region'] == region]
            
            for _, row in region_traits.iterrows():
                trait = row['Trait']
                
                # Add OR value
                region_data[f'OR_{trait}'] = row['OR']
                
                # Add confidence intervals
                region_data[f'CI_lower_{trait}'] = row['CI_lower']
                region_data[f'CI_upper_{trait}'] = row['CI_upper']
                
                # Add p-value
                region_data[f'p_value_{trait}'] = row['p_value']
            
            formatted_results.append(region_data)
        
        # Create formatted DataFrame
        formatted_df = pd.DataFrame(formatted_results)
        
        # Reorder columns to match image format: Region, then OR columns, then CI columns, then p-value columns
        # Get column order
        region_col = ['Region']
        or_cols = sorted([col for col in formatted_df.columns if col.startswith('OR_')])
        ci_lower_cols = sorted([col for col in formatted_df.columns if col.startswith('CI_lower_')])
        ci_upper_cols = sorted([col for col in formatted_df.columns if col.startswith('CI_upper_')])
        p_value_cols = sorted([col for col in formatted_df.columns if col.startswith('p_value_')])
        
        # Reorder columns
        formatted_df = formatted_df[region_col + or_cols + ci_lower_cols + ci_upper_cols + p_value_cols]
        
        # Save formatted results
        formatted_output = output_csv.replace('.csv', '_formatted.csv')
        formatted_df.to_csv(formatted_output, index=False, float_format='%.6f')
        print(f"Formatted results saved to: {formatted_output}")
        
        # Save detailed results
        detailed_output = output_csv.replace('.csv', '_detailed.csv')
        detailed_df.to_csv(detailed_output, index=False)
        print(f"Detailed results saved to: {detailed_output}")
        
        # Store results in class attribute
        self.results[analysis_type] = {
            'detailed': detailed_df,
            'formatted': formatted_df
        }
        
        return formatted_df, detailed_df
    
    def create_or_matrix(self, analysis_type=None, output_csv="OR_matrix.csv"):
        """
        Create OR matrix similar to the image format
        
        Parameters:
            analysis_type: Analysis type, if None returns the first analysis type
            output_csv: Path to save OR matrix
        """
        if not self.results:
            print("No analysis results found. Run batch_analyze_traits first.")
            return None
        
        if analysis_type is None:
            analysis_type = list(self.results.keys())[0]
        
        if analysis_type not in self.results:
            print(f"No results found for analysis type: {analysis_type}")
            return None
        
        formatted_df = self.results[analysis_type]['formatted']
        
        # Extract OR columns only
        or_cols = [col for col in formatted_df.columns if col.startswith('OR_')]
        or_matrix = formatted_df[['Region'] + or_cols].copy()
        
        # Save OR matrix
        or_matrix.to_csv(output_csv, index=False, float_format='%.6f')
        print(f"OR matrix saved to: {output_csv}")
        
        return or_matrix
    
    def create_ci_matrix(self, analysis_type=None, output_csv="CI_matrix.csv"):
        """
        Create confidence interval matrix
        
        Parameters:
            analysis_type: Analysis type, if None returns the first analysis type
            output_csv: Path to save CI matrix
        """
        if not self.results:
            print("No analysis results found. Run batch_analyze_traits first.")
            return None
        
        if analysis_type is None:
            analysis_type = list(self.results.keys())[0]
        
        if analysis_type not in self.results:
            print(f"No results found for analysis type: {analysis_type}")
            return None
        
        formatted_df = self.results[analysis_type]['formatted']
        
        # Extract CI columns
        ci_lower_cols = sorted([col for col in formatted_df.columns if col.startswith('CI_lower_')])
        ci_upper_cols = sorted([col for col in formatted_df.columns if col.startswith('CI_upper_')])
        
        ci_matrix = formatted_df[['Region'] + ci_lower_cols + ci_upper_cols].copy()
        
        # Save CI matrix
        ci_matrix.to_csv(output_csv, index=False, float_format='%.6e')
        print(f"CI matrix saved to: {output_csv}")
        
        return ci_matrix
    
    def plot_or_heatmap(self, analysis_type=None, figsize=(12, 8)):
        """
        Create a heatmap of OR values similar to the image
        
        Parameters:
            analysis_type: Analysis type, if None uses the first analysis type
            figsize: Figure size
        """
        if not self.results:
            print("No analysis results found. Run batch_analyze_traits first.")
            return
        
        if analysis_type is None:
            analysis_type = list(self.results.keys())[0]
        
        if analysis_type not in self.results:
            print(f"No results found for analysis type: {analysis_type}")
            return
        
        or_matrix = self.create_or_matrix(analysis_type)
        
        if or_matrix is None:
            return
        
        # Prepare data for heatmap
        heatmap_data = or_matrix.set_index('Region')
        
        # Create heatmap
        plt.figure(figsize=figsize)
        
        # Use a diverging colormap centered at 1 (no effect)
        cmap = plt.cm.RdBu_r
        center = 1
        
        # Create the heatmap
        im = plt.imshow(heatmap_data.values, cmap=cmap, aspect='auto', 
                       norm=plt.cm.colors.CenteredNorm(vcenter=center))
        
        # Add colorbar
        plt.colorbar(im, label='Odds Ratio (OR)')
        
        # Set ticks and labels
        plt.xticks(np.arange(len(heatmap_data.columns)), 
                  [col.replace('OR_', '') for col in heatmap_data.columns], 
                  rotation=45, ha='right')
        plt.yticks(np.arange(len(heatmap_data.index)), heatmap_data.index)
        
        # Add grid
        plt.grid(False)
        
        # Add title
        plt.title(f'OR Values Heatmap - {analysis_type}')
        
        # Add OR values as text
        for i in range(len(heatmap_data.index)):
            for j in range(len(heatmap_data.columns)):
                value = heatmap_data.iloc[i, j]
                if not pd.isna(value):
                    text = f'{value:.2f}' if value < 100 else f'{value:.0f}'
                    plt.text(j, i, text, ha='center', va='center', 
                            color='white' if abs(value - center) > 0.5 else 'black',
                            fontsize=8)
        
        plt.tight_layout()
        plt.show()
    
    def get_summary_statistics(self, analysis_type=None):
        """
        Get statistical summary of analysis results
        
        Parameters:
            analysis_type: Analysis type, if None returns results for all analyses
        """
        if analysis_type:
            if analysis_type not in self.results:
                print(f"No analysis results found for {analysis_type}")
                return None
            
            detailed_df = self.results[analysis_type]['detailed']
            summary = {
                'analysis_type': analysis_type,
                'total_tests': len(detailed_df),
                'significant_results': len(detailed_df[detailed_df['p_value'] < 0.05]),
                'mean_OR': detailed_df['OR'].mean(),
                'median_OR': detailed_df['OR'].median(),
                'min_OR': detailed_df['OR'].min(),
                'max_OR': detailed_df['OR'].max(),
                'OR_std': detailed_df['OR'].std()
            }
            return summary
        else:
            summary = {}
            for analysis_type, result_dict in self.results.items():
                detailed_df = result_dict['detailed']
                summary[analysis_type] = {
                    'total_tests': len(detailed_df),
                    'significant_results': len(detailed_df[detailed_df['p_value'] < 0.05]),
                    'mean_OR': detailed_df['OR'].mean(),
                    'median_OR': detailed_df['OR'].median()
                }
            return summary
    
    def save_results(self, output_dir="OR_results"):
        """
        Save all analysis results to specified directory
        
        Parameters:
            output_dir: Output directory path
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        for analysis_type, result_dict in self.results.items():
            # Save formatted results
            formatted_csv = output_dir / f"{analysis_type}_formatted.csv"
            result_dict['formatted'].to_csv(formatted_csv, index=False, float_format='%.6f')
            
            # Save detailed results
            detailed_csv = output_dir / f"{analysis_type}_detailed.csv"
            result_dict['detailed'].to_csv(detailed_csv, index=False)
            
            # Save OR matrix
            or_matrix_csv = output_dir / f"{analysis_type}_OR_matrix.csv"
            or_matrix = result_dict['formatted'][['Region'] + 
                        [col for col in result_dict['formatted'].columns if col.startswith('OR_')]]
            or_matrix.to_csv(or_matrix_csv, index=False, float_format='%.6f')
            
            # Save summary statistics
            summary = self.get_summary_statistics(analysis_type)
            if summary:
                summary_path = output_dir / f"{analysis_type}_summary.txt"
                with open(summary_path, 'w', encoding='utf-8') as f:
                    for key, value in summary.items():
                        f.write(f"{key}: {value}\n")
        
        print(f"All results saved to directory: {output_dir}")
        
        # Return paths to saved files
        saved_files = [str(p) for p in output_dir.iterdir() if p.is_file()]
        return saved_files


