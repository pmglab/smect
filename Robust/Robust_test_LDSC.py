import pandas as pd
import numpy as np
from scipy import stats
from typing import List, Dict, Union, Optional


def calculate_correlation_matrix(
    file_paths: List[str], 
    annotation_col: str = 'annotation', 
    p_value_col: str = 'p_cauchy',
    correlation_method: str = 'spearman',
    index_names: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Calculate correlation matrix of p-values from multiple CSV files.
    
    This function reads multiple CSV files, extracts p-value columns, and computes
    a correlation matrix between the p-values from different files.
    
    Args:
        file_paths: List of paths to CSV files containing p-value data.
        annotation_col: Name of the annotation column (default: 'annotation').
        p_value_col: Name of the p-value column (default: 'p_cauchy').
        correlation_method: Correlation method ('spearman' or 'pearson', default: 'spearman').
        index_names: Optional list of names for the rows/columns of the result matrix.
    
    Returns:
        DataFrame containing the correlation matrix.
    
    Raises:
        ValueError: If input parameters are invalid or files lack required columns.
        FileNotFoundError: If any specified file cannot be found.
        RuntimeError: If errors occur during file processing or data merging.
    
    Example:
        >>> file_paths = ["file1.csv", "file2.csv", "file3.csv", "file4.csv"]
        >>> corr_matrix = calculate_correlation_matrix(file_paths)
        >>> print(corr_matrix)
    """
    
    # Validate input parameters
    if not file_paths:
        raise ValueError("File path list cannot be empty")
    
    if correlation_method not in ['spearman', 'pearson']:
        raise ValueError("Correlation method must be 'spearman' or 'pearson'")
    
    # Read and preprocess data from all files
    data_frames = []
    
    for i, file_path in enumerate(file_paths):
        try:
            df = pd.read_csv(file_path)
            
            # Verify required columns exist
            if annotation_col not in df.columns or p_value_col not in df.columns:
                raise ValueError(f"File {file_path} lacks required columns")
            
            # Select and rename relevant columns
            df_p = df[[annotation_col, p_value_col]].rename(
                columns={p_value_col: f'result_df{i+1}_p'}
            )
            
            data_frames.append(df_p)
            
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        except Exception as e:
            raise RuntimeError(f"Error processing file {file_path}: {str(e)}")
    
    # Merge all DataFrames on annotation column
    try:
        combined_df = data_frames[0]
        for df in data_frames[1:]:
            combined_df = combined_df.merge(df, on=annotation_col, how='inner')
    except Exception as e:
        raise RuntimeError(f"Error merging DataFrames: {str(e)}")
    
    # Extract p-value columns for correlation calculation
    p_value_cols = [col for col in combined_df.columns if col.endswith('_p')]
    
    # Initialize correlation matrix
    n = len(p_value_cols)
    corr_matrix = pd.DataFrame(np.zeros((n, n)), 
                               columns=p_value_cols, 
                               index=p_value_cols)
    
    # Calculate pairwise correlations
    for i, col1 in enumerate(p_value_cols):
        for j, col2 in enumerate(p_value_cols):
            if i <= j:  # Avoid redundant calculations
                # Handle missing values
                valid_mask = (~combined_df[col1].isna()) & (~combined_df[col2].isna())
                valid_data = combined_df[valid_mask]
                
                if len(valid_data) > 1:  # Require at least 2 valid data points
                    if correlation_method == 'spearman':
                        corr, p_val = stats.spearmanr(valid_data[col1], valid_data[col2])
                    else:  # Pearson correlation
                        corr, p_val = stats.pearsonr(valid_data[col1], valid_data[col2])
                    
                    corr_matrix.loc[col1, col2] = corr
                    corr_matrix.loc[col2, col1] = corr
                else:
                    corr_matrix.loc[col1, col2] = np.nan
                    corr_matrix.loc[col2, col1] = np.nan
    
    # Apply custom index names if provided
    if index_names is None:
        index_names = [f'result_df{i+1}' for i in range(n)]
    elif len(index_names) != n:
        raise ValueError(f"index_names length must be {n}")
    
    corr_matrix.index = index_names
    corr_matrix.columns = index_names
    
    return corr_matrix


def correlation_analysis(
    file_paths: List[str],
    **kwargs
) -> Dict[str, Union[pd.DataFrame, pd.Series]]:
    """
    Perform comprehensive correlation analysis on multiple data files.
    
    This function provides a complete analysis workflow including correlation matrix
    calculation and descriptive statistics for the combined dataset.
    
    Args:
        file_paths: List of paths to CSV files containing p-value data.
        **kwargs: Additional arguments passed to calculate_correlation_matrix.
    
    Returns:
        Dictionary containing:
            - correlation_matrix: DataFrame with correlation coefficients
            - descriptive_stats: Summary statistics for p-values
            - combined_data: Merged DataFrame containing all p-values
    
    Example:
        >>> results = correlation_analysis(file_paths, correlation_method='spearman')
        >>> print(results['correlation_matrix'])
        >>> print(results['descriptive_stats'])
    """
    
    # Calculate correlation matrix
    corr_matrix = calculate_correlation_matrix(file_paths, **kwargs)
    
    # Calculate descriptive statistics
    data_frames = []
    for i, file_path in enumerate(file_paths):
        df = pd.read_csv(file_path)
        p_col = kwargs.get('p_value_col', 'p_cauchy')
        if p_col in df.columns:
            df_p = df[[kwargs.get('annotation_col', 'annotation'), p_col]].rename(
                columns={p_col: f'result_df{i+1}_p'}
            )
            data_frames.append(df_p)
    
    # Merge data for descriptive statistics
    combined_df = data_frames[0]
    for df in data_frames[1:]:
        combined_df = combined_df.merge(df, on=kwargs.get('annotation_col', 'annotation'), how='inner')
    
    p_value_cols = [col for col in combined_df.columns if col.endswith('_p')]
    descriptive_stats = combined_df[p_value_cols].describe()
    
    return {
        'correlation_matrix': corr_matrix,
        'descriptive_stats': descriptive_stats,
        'combined_data': combined_df
    }


