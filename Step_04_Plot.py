import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scienceplots
import numpy as np
import os
from pathlib import Path
import glob

import argparse

# Set up the plotting style for better-looking figures
# plt.style.use('seaborn-v0_8')
# sns.set_palette("husl")
plt.style.use('science')

# Define the methods and their colors for consistent visualization
methods = ['CASPO', 'VNS', 'GA', 'ILP']
colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00']  # Red, Blue, Green, Orange (ColorBrewer Set1)
method_colors = dict(zip(methods, colors))
    
# Helper to plot a metric
def plot_metric(df, ax, metric, title, marker):

    ax.set_xlim(.0, 1.)
    ax.set_ylim(df[metric].min()-0.05, df[metric].max() + 0.05)  
    for m in methods:
        mdata = df[df['method'] == m]
        ax.plot(
            mdata['change_percent'], mdata[metric],
            marker=marker, linewidth=1.5, markersize=3,
            label=m, color=method_colors[m]
        )
    ax.set_xlabel('Change Percentage')
    ax.set_ylabel(title)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
        
        
def create_comparison_plots(df, dataset_name='toy'):
    """
    Create comprehensive comparison plots for different methods across change percentages.
    
    Parameters:
    -----------
    dataset_name : str
        Name of the dataset (used for folder creation)
    """
    
    # Create output directory
    output_dir = Path(f"output/{dataset_name}")
    output_dir.mkdir(exist_ok=True)
    
    # Define key metrics to focus on
    key_metrics = {
        'Similarity Metrics': ['jaccard_similarity', 'hamming_similarity', 'composite_score'],
        'Coverage Metrics': ['precision', 'recall', 'f1_score'],
    }
    
    # Create figure 1: Core similarity metrics comparison
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    fig.suptitle(f'Core Similarity Metrics Comparison Across Methods of {dataset_name}', fontsize=16, fontweight='bold')

    plot_metric(df, axes[0, 0], 'jaccard_similarity', 'Jaccard Similarity Performance', 'o')
    plot_metric(df, axes[0, 1], 'hamming_similarity', 'Hamming Similarity Performance', 's')
    plot_metric(df, axes[1, 0], 'jaccard_topology', 'Jaccard Topology Performance', '^')
    plot_metric(df, axes[1, 1], 'f1_score', 'F1 Score Performance', 'D')
    plt.tight_layout()
    fig.savefig(output_dir / 'core_similarity_metrics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 2: Coverage and precision metrics
    fig, axes = plt.subplots(2, 1, figsize=(8, 6))
    fig.suptitle(f'Coverage and Precision Metrics Analysis of {dataset_name}', fontsize=16, fontweight='bold')

    plot_metric(df, axes[0], 'precision', 'Precision Performance', 'o')
    plot_metric(df, axes[1], 'recall', 'Recall Performance', 's')
    plt.tight_layout()
    fig.savefig(output_dir / 'coverage_precision_metrics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
            
    # Create figure 3: Performance summary heatmap
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Calculate average performance for each method across all change percentages
    performance_summary = df.groupby('method')[['jaccard_similarity', 'hamming_similarity', 'jaccard_topology',
                                                'composite_score', 'f1_score', 'precision', 'recall']].mean()

    sns.heatmap(
        performance_summary.T, annot=True, fmt='.3f',
        cmap='RdYlBu_r', center=0.5, annot_kws={'color': 'black', 'fontsize': 9},  
        cbar_kws={'label': 'Performance Score'},  
        linewidths=0.5,  linecolor='gray',   
        ax=ax
    )
    ax.set_title(f'Average Performance Summary Across All Metrics of {dataset_name}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Methods')
    ax.set_ylabel('Metrics')
    plt.tight_layout()
    fig.savefig(output_dir / 'performance_summary_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 4: jaccard_topology plot    
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_metric(df, ax, 'jaccard_topology', 'Jaccard Topology Performance', 'o')
    plt.tight_layout()
    fig.savefig(output_dir / 'jaccard_topology.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 5: Runtime comparison
    if 'total_time' in df.columns:      
        fig, ax = plt.subplots(figsize=(8, 6))  
        
        for method in methods:
            method_data = df[df['method'] == method]
            ax.semilogy(method_data['change_percent'], method_data['total_time'], 
                       marker='o', linewidth=1.5, markersize=3, label=method, color=method_colors[method])
        
        ax.set_xlabel('Change Percentage')
        ax.set_ylabel('Total Time (log scale)')
        ax.set_title(f'Runtime Comparison Across Methods of {dataset_name}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'runtime_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Create figure 6: Robustness analysis (performance degradation)    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    # Calculate performance degradation (relative to 0.0 change percentage)
    metrics_to_analyze = ['composite_score', 'f1_score', 'jaccard_similarity', 'hamming_similarity', 'jaccard_topology']

    for i, method in enumerate(methods):
        ax_sub = axes[i // 2, i % 2]
        ax_sub.set_xlim(0., 1.)
        ax_sub.set_ylim(0.0, 1.5)
        method_data = df[df['method'] == method]
        for j, metric in enumerate(metrics_to_analyze):
            
            baseline = method_data[method_data['change_percent'] == 0.0][metric].iloc[0]
            
            # Calculate relative performance (performance / baseline)
            relative_performance = method_data[metric] / baseline if baseline != 0 else method_data[metric]
            
            ax_sub.plot(method_data['change_percent'], relative_performance, 
                       marker='o', linewidth=1.5, markersize=3, label=metrics_to_analyze[j].replace("_", " "))

        ax_sub.set_xlabel('Change Percentage')
        ax_sub.set_ylabel(f'Relative {method.replace("_", " ").upper()}')
        ax_sub.set_title(f'Robustness: {method.replace("_", " ").upper()} of {dataset_name}')
        ax_sub.legend()
        ax_sub.grid(True, alpha=0.3)
        ax_sub.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'robustness_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 7 plot number of reconstruct attractors
    fig, ax = plt.subplots(figsize=(8, 6))

    for method in methods:
        method_data = df[df['method'] == method]
        ax.semilogy(method_data['change_percent'], method_data['recon_total'], 
                    marker='o', linewidth=1.5, markersize=3, label=method, color=method_colors[method])

    ax.set_xlabel('Change Percentage')
    ax.set_ylabel('Number of Reconstructed Attractors')
    ax.set_title(f'Number of Reconstructed Attractors Across Methods of {dataset_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'reconstructed_attractors.png', dpi=300, bbox_inches='tight')
    plt.close()            
    
    print(f"All plots have been saved to the '{dataset_name}' directory!")
    print(f"Generated plots:")
    print(f"  - core_similarity_metrics.png")
    print(f"  - coverage_precision_metrics.png") 
    print(f"  - performance_summary_heatmap.png")
    print(f"  - jaccard_topology.png")
    if 'total_time' in df.columns:
        print(f"  - runtime_comparison.png")
    print(f"  - robustness_analysis.png")
    print(f"  - reconstructed_attractors.png")

def compute_attractor_count_penalty(self, expected_attractors=None):
    """
    Penalize networks that deviate significantly from expected attractor counts.
    """
    if expected_attractors is None:
        expected_attractors = len(self.original_attractors)
    
    actual_count = len(self.reconstructed_attractors)
    
    # Use a sigmoid-like penalty that's gentle for small deviations
    # but harsh for large ones
    ratio = actual_count / max(expected_attractors, 1)
    
    if ratio <= 1:
        penalty = 1.0  # No penalty for fewer attractors
    else:
        # Penalty increases rapidly after 2x the expected count
        penalty = 1.0 / (1.0 + 0.5 * (ratio - 1)**2)
    
    return penalty

def rescale_similarity(df, metric):
    recon_size = df['recon_total']
    orig_size = df['orig_total']
    ratio = orig_size / max(orig_size, 1)
    if ratio <= 1:
        penalty = 1.0  # No penalty for fewer attractors
    else:
        # Penalty increases rapidly after 2x the expected count
        penalty = 2.0 / (1.0 + np.exp(0.5 * (ratio - 2)))
    return df[metric] * penalty

def load_results(pattern):
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No files match pattern {pattern}")
    dfs = []
    for p in paths:
        df = pd.read_csv(p)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def load_and_plot_results(dataset_name='toy'):
    
    pattern = f"output/comparison_{dataset_name}_*.csv"  # note: `*` matches any characters except `/`

    full_df = load_results(pattern)
    print(full_df.head())
    full_df['jaccard_similarity'] = full_df.apply(lambda row: rescale_similarity(row, 'jaccard_similarity'), axis=1)
    full_df['hamming_similarity'] = full_df.apply(lambda row: rescale_similarity(row, 'hamming_similarity'), axis=1)    
    
    print(full_df.head())
    
    avg_columns = [
        'jaccard_similarity', 'hamming_similarity', 'composite_score',
        'precision', 'recall', 'f1_score', 'true_positives', 'orig_total', 'recon_total',
        'total_time', 'jaccard_topology'
    ]

    group_columns = ['method', 'change_percent']

    df = full_df.groupby(group_columns, dropna=False)[avg_columns].mean().reset_index()
    
    df['dataset'] = dataset_name
    # Create all plots
    create_comparison_plots(df, dataset_name)
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print("=" * 50)
    for method in df['method'].unique():
        method_data = df[df['method'] == method]
        avg_topo    = method_data['jaccard_topology'].mean()
        avg_composite = method_data['composite_score'].mean()
        avg_jaccard = method_data['jaccard_similarity'].mean()
        avg_hamming = method_data['hamming_similarity'].mean()
        avg_f1      = method_data['f1_score'].mean()
        print(f"{method.upper():>6}: Avg Composite Score = {avg_composite:.3f}, Avg Jaccard = {avg_jaccard:.3f}, Avg Hamming = {avg_hamming:.3f}, Avg F1 score = {avg_f1:.3f}, Avg Topology = {avg_topo:.3f}")

    toy_file = glob.glob(f"output/comparison_toy_*.csv")
    if dataset_name != 'toy' and (len(toy_file) > 0): 
        print("Loading toy dataset for comparison...")       
        toy_df = load_results(f"output/comparison_toy_*.csv")
        toy_df = toy_df.groupby(group_columns, dropna=False)[avg_columns].mean().reset_index()
        toy_df['dataset'] = "Toy"
        combined_df = pd.concat([df, toy_df], ignore_index=True)

        for metric in ['jaccard_similarity', 'hamming_similarity', 'composite_score', 'total_time', 'jaccard_topology', 'recon_total']:
            fig, axes = plt.subplots(2, 2, figsize=(8, 6))
            fig.suptitle(f'{metric} Similarity Metrics Comparison Across Network', fontsize=16, fontweight='bold')
            ymin, ymax = combined_df[metric].min(), combined_df[metric].max()
                    

            for j, dataset in enumerate([dataset_name, 'Toy']):
                dataset_data = combined_df[combined_df['dataset'] == dataset]
                for i, m in enumerate(methods):
                    ax = axes[i%2, i//2]
                    ax.set_xlim(.0, 1.)
                    ax.set_ylim(ymin - 0.05, ymax + 0.05)
                    method_data = dataset_data[dataset_data['method'] == m]
                    ax.plot(
                        method_data['change_percent'], 
                        method_data[metric],
                        marker='o', 
                        linewidth=1.5,
                        markersize=3,
                        label=dataset,
                    )
                    ax.set_xlabel('Change Percentage')
                    ax.set_ylabel(f"{metric} Comparison")
                    ax.set_title(f'Method: {m.upper()}')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
            plt.tight_layout()
            fig.savefig(f'output/{dataset_name}/size_comparison_{metric}.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

# Example usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run model comparisons multiple times on a given dataset"
    )
    parser.add_argument(
        "-d", "--dataset",
        type=str,
        default="toy",
        help="Dataset to use for comparisons (default: 'toy')"
    )
    args = parser.parse_args()
    dataset = args.dataset
    # Load and create plots
    load_and_plot_results(dataset_name=dataset)
    