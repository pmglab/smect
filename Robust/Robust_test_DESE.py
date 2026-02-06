import numpy as np
import pandas as pd
from scipy.stats import cauchy, spearmanr


def cauchy_combination_dese_test(pvals, weights=None):
    """
    Calculate Cauchy combination test statistic and combined p-value

    Parameters:
        pvals (list/array): List of p-values
        weights (list/array): List of weights (default: equal weights)

    Returns:
        tuple: (Test statistic T, Combined p-value)
    """
    pvals = np.asarray(pvals)
    if np.any(pvals < 0) or np.any(pvals > 1):
        raise ValueError("P-values must be in the range [0, 1]")

    n = len(pvals)
    weights = np.ones(n) / n if weights is None else np.asarray(weights)

    # Transform p-values: Handle boundary values (0 or 1) to avoid infinity
    safe_pvals = np.clip(pvals, 1e-15, 1 - 1e-15)
    transformed = np.tan((0.5 - safe_pvals) * np.pi)

    T = np.sum(weights * transformed)  # Test statistic
    p_combined = cauchy.sf(T)  # Calculate p-value using Cauchy survival function

    return T, p_combined


def process_enrichment_data(enrichment_paths, truth_paths, layers=None):
    """
    Process enrichment analysis data and calculate Cauchy combination test

    Parameters:
        enrichment_paths (list): List of enrichment analysis file paths
        truth_paths (list): List of ground truth annotation file paths
        layers (list): List of brain region layers, default: 7 standard layers

    Returns:
        dict: Dictionary containing result DataFrames for each sample
    """
    if layers is None:
        layers = ['Layer_1', 'Layer_2', 'Layer_3', 'Layer_4', 'Layer_5', 'Layer_6', 'WM']

    if len(enrichment_paths) != len(truth_paths):
        raise ValueError("Number of enrichment analysis file paths and ground truth file paths must be equal")

    results = {}

    for i, (enrich_path, truth_path) in enumerate(zip(enrichment_paths, truth_paths)):
        # Read data
        df = pd.read_excel(enrich_path, engine="xlrd")
        ann = pd.read_csv(truth_path, delimiter='\t')

        # Data preprocessing
        df = df.drop(df.index[-2:])
        ann.columns = ['Condition', 'annotation']
        ann = ann.set_index('Condition')
        df = df.set_index('Condition')

        # Merge data
        df = pd.concat([df, ann], axis=1)
        df = df.dropna(subset=['annotation'])

        # Select relevant columns and calculate p-values
        dfp = df.iloc[:, [1, 3]]
        dfp['p_value'] = dfp['EnrichmentScore']

        # Calculate Cauchy combination test for each brain region layer
        layer_results = []
        for layer in layers:
            # Extract p-values for current brain region layer (non-NaN)
            mask = (dfp['annotation'] == layer) & (~dfp['EnrichmentScore'].isna())
            layer_pvals = dfp.loc[mask, 'p_value'].tolist()

            if layer_pvals:
                T, p_combined = cauchy_combination_dese_test(layer_pvals)
                num_pvals = len(layer_pvals)
            else:
                T, p_combined, num_pvals = np.nan, np.nan, 0

            layer_results.append({
                'Annotation': layer,
                'Cauchy_Stat': T,
                'Combined_p_value': p_combined,
                'Num_p_values': num_pvals
            })

        results[f'df{i + 1}'] = pd.DataFrame(layer_results)

    return results


def calculate_correlation_matrix(results_dict):
    """
    Calculate correlation coefficient matrix for negative log-transformed p-values across multiple result DataFrames

    Parameters:
        results_dict (dict): Dictionary containing result DataFrames

    Returns:
        pd.DataFrame: Correlation coefficient matrix of negative log p-values
    """
    # Merge all results
    combined_df = None
    for i, (key, df) in enumerate(results_dict.items()):
        temp_df = df[['Annotation', 'Combined_p_value']].copy()
        temp_df = temp_df.rename(columns={'Combined_p_value': f'df{i + 1}_p'})

        if combined_df is None:
            combined_df = temp_df
        else:
            combined_df = combined_df.merge(temp_df, on='Annotation')

    # Extract p-value columns and convert to negative log10(p-values)
    p_value_cols = [col for col in combined_df.columns if '_p' in col]

    # Convert p-values to negative log10 scale
    for col in p_value_cols:
        # Avoid log(0) by clipping very small p-values and handle NaN values
        combined_df[f'{col}_neg_log10'] = -np.log10(np.clip(combined_df[col], 1e-16, 1.0))

    # Get the negative log transformed columns
    neg_log_cols = [col for col in combined_df.columns if 'neg_log10' in col]

    n = len(neg_log_cols)
    corr_matrix = pd.DataFrame(np.zeros((n, n)),
                               columns=neg_log_cols,
                               index=neg_log_cols)

    # Calculate Spearman correlation for negative log transformed p-values
    for i, col1 in enumerate(neg_log_cols):
        for j, col2 in enumerate(neg_log_cols):
            if i <= j:
                valid_mask = (~combined_df[col1].isna()) & (~combined_df[col2].isna())
                if valid_mask.sum() > 1:  # At least 2 valid values required
                    corr, p_val = spearmanr(combined_df.loc[valid_mask, col1],
                                            combined_df.loc[valid_mask, col2])
                    corr_matrix.loc[col1, col2] = corr
                    corr_matrix.loc[col2, col1] = corr
                else:
                    corr_matrix.loc[col1, col2] = np.nan
                    corr_matrix.loc[col2, col1] = np.nan

    # Simplify output formatting
    simplified_index = [f'df{i + 1}_neg_log_p' for i in range(n)]
    corr_matrix.index = simplified_index
    corr_matrix.columns = simplified_index

    return corr_matrix


def main_analysis():
    """
    Main analysis function: Execute complete analysis pipeline
    """
    # Define file paths
    enrichment_paths = [
        r"C:/Users/pmgl903/Desktop/kgg/kggsum3/SCZ/kggsum/HS_151507_GSS.txt.gz.enrichment.xls",
        r"C:/Users/pmgl903/Desktop/kgg/kggsum3/SCZ/kggsum/HS_151508_GSS.txt.gz.enrichment.xls",
        r"C:/Users/pmgl903/Desktop/kgg/kggsum3/SCZ/kggsum/HS_151509_GSS.txt.gz.enrichment.xls",
        r"C:/Users/pmgl903/Desktop/kgg/kggsum3/SCZ/kggsum/HS_151510_GSS.txt.gz.enrichment.xls"
    ]

    truth_paths = [
        'C:/Users/pmgl903/Desktop/stprojct/st/MNMST-Data/MNMST-Data/DLPFC/151507/151507_truth.txt',
        'C:/Users/pmgl903/Desktop/stprojct/st/MNMST-Data/MNMST-Data/DLPFC/151508/151508_truth.txt',
        'C:/Users/pmgl903/Desktop/stprojct/st/MNMST-Data/MNMST-Data/DLPFC/151509/151509_truth.txt',
        'C:/Users/pmgl903/Desktop/stprojct/st/MNMST-Data/MNMST-Data/DLPFC/151510/151510_truth.txt'
    ]

    # Process data and calculate Cauchy combination test
    results = process_enrichment_data(enrichment_paths, truth_paths)

    # Calculate correlation matrix using negative log transformed p-values
    corr_matrix = calculate_correlation_matrix(results)

    # Print results
    print("Cauchy combination test results:")
    for key, df in results.items():
        print(f"\n{key}:")
        print(df)

    print("\nPairwise Spearman correlation matrix for Negative Log10 Transformed Combined_p_value:")
    print(corr_matrix)

    return results, corr_matrix


# Execute the main analysis if run as script
if __name__ == "__main__":
    results, corr_matrix = main_analysis()
