import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path
import glob

import argparse

# Set up the plotting style for better-looking figures
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# Define the methods and their colors for consistent visualization
methods = ['Caspo', 'VNS', 'GA', 'ILP']
colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00']  # Red, Blue, Green, Orange (ColorBrewer Set1)
method_colors = dict(zip(methods, colors))
    
    
# Helper to plot a metric
def plot_metric(df, ax, metric, title, marker):
    for m in methods:
        mdata = df[df['method'] == m]
        ax.plot(
            mdata['change_percent'], mdata[metric],
            marker=marker, linewidth=2.5, markersize=6,
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
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Core Similarity Metrics Comparison Across Methods', fontsize=16, fontweight='bold')

    plot_metric(df, axes[0, 0], 'jaccard_similarity', 'Jaccard Similarity Performance', 'o')
    plot_metric(df, axes[0, 1], 'hamming_similarity', 'Hamming Similarity Performance', 's')
    plot_metric(df, axes[1, 0], 'composite_score', 'Overall Composite Score', '^')
    plot_metric(df, axes[1, 1], 'f1_score', 'F1 Score Performance', 'D')
    plt.tight_layout()
    fig.savefig(output_dir / 'core_similarity_metrics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 2: Coverage and precision metrics
    fig, axes = plt.subplots(2, 1, figsize=(15, 12))
    fig.suptitle('Coverage and Precision Metrics Analysis', fontsize=16, fontweight='bold')

    plot_metric(df, axes[0], 'precision', 'Precision Performance', 'o')
    plot_metric(df, axes[1], 'recall', 'Recall Performance', 's')
    plt.tight_layout()
    fig.savefig(output_dir / 'coverage_precision_metrics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
            
    # Create figure 3: Performance summary heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Calculate average performance for each method across all change percentages
    performance_summary = df.groupby('method')[['jaccard_similarity', 'hamming_similarity', 'jaccard_topology',
                                                'composite_score', 'f1_score', 'precision', 'recall']].mean()

    sns.heatmap(
        performance_summary.T, annot=True, fmt='.3f',
        cmap='RdYlBu_r', center=0.5,
        cbar_kws={'label': 'Performance Score'},
        ax=ax
    )
    ax.set_title('Average Performance Summary Across All Metrics', fontsize=14, fontweight='bold')
    ax.set_xlabel('Methods')
    ax.set_ylabel('Metrics')
    plt.tight_layout()
    fig.savefig(output_dir / 'performance_summary_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 4: jaccard_topology plot    
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_metric(df, ax, 'jaccard_topology', 'Jaccard Topology Performance', 'o')
    plt.tight_layout()
    fig.savefig(output_dir / 'jaccard_topology.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Create figure 5: Runtime comparison
    if 'total_time' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for method in methods:
            method_data = df[df['method'] == method]
            ax.semilogy(method_data['change_percent'], method_data['total_time'], 
                       marker='o', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
        
        ax.set_xlabel('Change Percentage')
        ax.set_ylabel('Total Time (log scale)')
        ax.set_title('Runtime Comparison Across Methods')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'runtime_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Create figure 6: Robustness analysis (performance degradation)
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Calculate performance degradation (relative to 0.0 change percentage)
    metrics_to_analyze = ['composite_score', 'jaccard_similarity', 'f1_score']
    
    for i, metric in enumerate(metrics_to_analyze):
        ax_sub = plt.subplot(2, 2, i+1)
        
        for method in methods:
            method_data = df[df['method'] == method].sort_values('change_percent')
            baseline = method_data[method_data['change_percent'] == 0.0][metric].iloc[0]
            
            # Calculate relative performance (performance / baseline)
            relative_performance = method_data[metric] / baseline if baseline != 0 else method_data[metric]
            
            ax_sub.plot(method_data['change_percent'], relative_performance, 
                       marker='o', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
        
        ax_sub.set_xlabel('Change Percentage')
        ax_sub.set_ylabel(f'Relative {metric.replace("_", " ").title()}')
        ax_sub.set_title(f'Robustness: {metric.replace("_", " ").title()}')
        ax_sub.legend()
        ax_sub.grid(True, alpha=0.3)
        ax_sub.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'robustness_analysis.png', dpi=300, bbox_inches='tight')
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

def load_and_plot_results(dataset_name='toy'):
    
    pattern = f"output/comparison_{dataset_name}_*.csv"  # note: `*` matches any characters except `/`

    matching_files = glob.glob(pattern)   
    
    # Read and concatenate all CSV files
    df_list = []
    for path in matching_files:
        df = pd.read_csv(path)
        df_list.append(df)

    full_df = pd.concat(df_list, ignore_index=True)
    
    avg_columns = [
        'jaccard_similarity', 'hamming_similarity', 'composite_score',
        'precision', 'recall', 'f1_score', 'true_positives', 'orig_total', 'recon_total',
        'total_time', 'jaccard_topology'
    ]

    group_columns = ['method', 'change_percent']

    df = full_df.groupby(group_columns, dropna=False)[avg_columns].mean().reset_index()

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
    