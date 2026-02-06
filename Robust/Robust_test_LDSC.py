import pandas as pd
import numpy as np
from scipy import stats
from typing import List, Dict, Union, Optional


def calculate_LDSC_correlation_matrix(
        file_paths: List[str],
        annotation_col: str = 'annotation',
        p_value_col: str = 'p_cauchy',
        correlation_method: str = 'spearman',
        index_names: Optional[List[str]] = None,
        apply_neg_log_transform: bool = True  # New parameter: whether to apply negative log transformation
) -> pd.DataFrame:
    """
    Calculate correlation matrix of p-values from multiple CSV files

    Parameters:
        file_paths (List[str]): List of CSV file paths
        annotation_col (str): Annotation column name, defaults to 'annotation'
        p_value_col (str): P-value column name, defaults to 'p_cauchy'
        correlation_method (str): Correlation calculation method, 'spearman' or 'pearson', defaults to 'spearman'
        index_names (Optional[List[str]]): Row/column names for result matrix, defaults to None
        apply_neg_log_transform (bool): Whether to apply negative log transformation to p-values, defaults to True

    Returns:
        pd.DataFrame: Correlation matrix

    Example:
        file_paths = [
            "path/to/file1.csv",
            "path/to/file2.csv",
            "path/to/file3.csv",
            "path/to/file4.csv"
        ]
        corr_matrix = calculate_correlation_matrix(file_paths)
    """

    # Validate input parameters
    if not file_paths:
        raise ValueError("File path list cannot be empty")

    if correlation_method not in ['spearman', 'pearson']:
        raise ValueError("Correlation method must be 'spearman' or 'pearson'")

    # Read and preprocess data
    data_frames = []

    for i, file_path in enumerate(file_paths):
        try:
            # Read CSV file
            df = pd.read_csv(file_path)

            # Validate required columns exist
            if annotation_col not in df.columns or p_value_col not in df.columns:
                raise ValueError(f"File {file_path} is missing required columns")

            # Apply negative log transformation to p-values (if enabled)
            if apply_neg_log_transform:
                # Handle p-values of 0 to avoid log calculation errors[6](@ref)
                p_values = df[p_value_col].copy()
                # Replace 0 values with very small positive numbers (close to machine precision)
                p_values[p_values == 0] = np.finfo(float).eps
                # Calculate negative log: -log10(p) or -ln(p), using natural log here[6](@ref)
                transformed_p = -np.log(p_values)
                transformed_col_name = f'result_df{i + 1}_neg_log_p'

                # Create new DataFrame with transformed p-values
                df_p = pd.DataFrame({
                    annotation_col: df[annotation_col],
                    transformed_col_name: transformed_p
                })
            else:
                # No transformation, use original p-values
                df_p = df[[annotation_col, p_value_col]].rename(
                    columns={p_value_col: f'result_df{i + 1}_p'}
                )

            data_frames.append(df_p)

        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        except Exception as e:
            raise RuntimeError(f"Error processing file {file_path}: {str(e)}")

    # Merge all DataFrames
    try:
        combined_df = data_frames[0]
        for df in data_frames[1:]:
            combined_df = combined_df.merge(df, on=annotation_col, how='inner')
    except Exception as e:
        raise RuntimeError(f"Error merging DataFrames: {str(e)}")

    # Extract p-value columns (select corresponding column names based on transformation)
    if apply_neg_log_transform:
        p_value_cols = [col for col in combined_df.columns if col.endswith('_neg_log_p')]
    else:
        p_value_cols = [col for col in combined_df.columns if col.endswith('_p')]

    # Calculate correlation matrix[1,2](@ref)
    n = len(p_value_cols)
    corr_matrix = pd.DataFrame(np.zeros((n, n)),
                               columns=p_value_cols,
                               index=p_value_cols)

    for i, col1 in enumerate(p_value_cols):
        for j, col2 in enumerate(p_value_cols):
            if i <= j:  # Avoid duplicate calculations
                # Handle missing values
                valid_mask = (~combined_df[col1].isna()) & (~combined_df[col2].isna())
                valid_data = combined_df[valid_mask]

                if len(valid_data) > 1:  # Need at least 2 valid data points
                    if correlation_method == 'spearman':
                        corr, p_val = stats.spearmanr(valid_data[col1], valid_data[col2])
                    else:  # pearson
                        corr, p_val = stats.pearsonr(valid_data[col1], valid_data[col2])

                    corr_matrix.loc[col1, col2] = corr
                    corr_matrix.loc[col2, col1] = corr
                else:
                    corr_matrix.loc[col1, col2] = np.nan
                    corr_matrix.loc[col2, col1] = np.nan

    # Beautify output
    if index_names is None:
        index_names = [f'result_df{i + 1}' for i in range(n)]
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
    Perform comprehensive correlation analysis, returning correlation matrix and descriptive statistics

    Parameters:
        file_paths (List[str]): List of file paths
        **kwargs: Parameters passed to calculate_correlation_matrix

    Returns:
        Dict: Dictionary containing correlation matrix and descriptive statistics
    """
    # Calculate correlation matrix
    corr_matrix = calculate_correlation_matrix(file_paths, **kwargs)

    # Calculate descriptive statistics
    data_frames = []
    for i, file_path in enumerate(file_paths):
        df = pd.read_csv(file_path)
        p_col = kwargs.get('p_value_col', 'p_cauchy')
        apply_transform = kwargs.get('apply_neg_log_transform', True)

        if p_col in df.columns:
            if apply_transform:
                p_values = df[p_col].copy()
                p_values[p_values == 0] = np.finfo(float).eps
                transformed_p = -np.log(p_values)
                transformed_col_name = f'result_df{i + 1}_neg_log_p'

                df_p = pd.DataFrame({
                    kwargs.get('annotation_col', 'annotation'): df[kwargs.get('annotation_col', 'annotation')],
                    transformed_col_name: transformed_p
                })
            else:
                df_p = df[[kwargs.get('annotation_col', 'annotation'), p_col]].rename(
                    columns={p_col: f'result_df{i + 1}_p'}
                )
            data_frames.append(df_p)

    combined_df = data_frames[0]
    for df in data_frames[1:]:
        combined_df = combined_df.merge(df, on=kwargs.get('annotation_col', 'annotation'), how='inner')

    if apply_transform:
        p_value_cols = [col for col in combined_df.columns if col.endswith('_neg_log_p')]
    else:
        p_value_cols = [col for col in combined_df.columns if col.endswith('_p')]

    descriptive_stats = combined_df[p_value_cols].describe()

    return {
        'correlation_matrix': corr_matrix,
        'descriptive_stats': descriptive_stats,
        'combined_data': combined_df
    }
