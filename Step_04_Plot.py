import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path

# Set up the plotting style for better-looking figures
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def create_comparison_plots(df, dataset_name='toy'):
    """
    Create comprehensive comparison plots for different methods across change percentages.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the results with all the specified columns
    dataset_name : str
        Name of the dataset (used for folder creation)
    """
    
    # Create output directory
    output_dir = Path(f"output/{dataset_name}")
    output_dir.mkdir(exist_ok=True)
    
    # Define the methods and their colors for consistent visualization
    methods = ['Caspo', 'VNS', 'GA', 'ILP']
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
    method_colors = dict(zip(methods, colors))
    
    # Define key metrics to focus on
    key_metrics = {
        'Similarity Metrics': ['jaccard_similarity', 'hamming_similarity', 'composite_score'],
        'Coverage Metrics': ['precision', 'recall', 'f1_score'],
        'Stability Metrics': ['stability_correlation', 'size_ratio'],
        'Functional Metrics': ['functional_similarity', 'pattern_overlap']
    }
    
    # Create figure 1: Core similarity metrics comparison
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Core Similarity Metrics Comparison Across Methods', fontsize=16, fontweight='bold')
    
    # Plot Jaccard Similarity
    ax1 = axes[0, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax1.plot(method_data['change_percent'], method_data['jaccard_similarity'], 
                marker='o', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax1.set_xlabel('Change Percentage')
    ax1.set_ylabel('Jaccard Similarity')
    ax1.set_title('Jaccard Similarity Performance')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot Hamming Similarity
    ax2 = axes[0, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax2.plot(method_data['change_percent'], method_data['hamming_similarity'], 
                marker='s', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax2.set_xlabel('Change Percentage')
    ax2.set_ylabel('Hamming Similarity')
    ax2.set_title('Hamming Similarity Performance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot Composite Score
    ax3 = axes[1, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax3.plot(method_data['change_percent'], method_data['composite_score'], 
                marker='^', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax3.set_xlabel('Change Percentage')
    ax3.set_ylabel('Composite Score')
    ax3.set_title('Overall Composite Score')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot F1 Score
    ax4 = axes[1, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax4.plot(method_data['change_percent'], method_data['f1_score'], 
                marker='D', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax4.set_xlabel('Change Percentage')
    ax4.set_ylabel('F1 Score')
    ax4.set_title('F1 Score Performance')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'core_similarity_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure 2: Coverage and precision metrics
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Coverage and Precision Metrics Analysis', fontsize=16, fontweight='bold')
    
    # Plot Precision
    ax1 = axes[0, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax1.plot(method_data['change_percent'], method_data['precision'], 
                marker='o', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax1.set_xlabel('Change Percentage')
    ax1.set_ylabel('Precision')
    ax1.set_title('Precision Performance')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot Recall
    ax2 = axes[0, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax2.plot(method_data['change_percent'], method_data['recall'], 
                marker='s', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax2.set_xlabel('Change Percentage')
    ax2.set_ylabel('Recall')
    ax2.set_title('Recall Performance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot Exact Matches
    ax3 = axes[1, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax3.plot(method_data['change_percent'], method_data['exact_matches'], 
                marker='^', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax3.set_xlabel('Change Percentage')
    ax3.set_ylabel('Exact Matches')
    ax3.set_title('Number of Exact Matches')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot Common Nodes
    ax4 = axes[1, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax4.plot(method_data['change_percent'], method_data['common_nodes'], 
                marker='D', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax4.set_xlabel('Change Percentage')
    ax4.set_ylabel('Common Nodes')
    ax4.set_title('Number of Common Nodes')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'coverage_precision_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure 3: Stability and functional metrics
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Stability and Functional Metrics Analysis', fontsize=16, fontweight='bold')
    
    # Plot Stability Correlation
    ax1 = axes[0, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax1.plot(method_data['change_percent'], method_data['stability_correlation'], 
                marker='o', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax1.set_xlabel('Change Percentage')
    ax1.set_ylabel('Stability Correlation')
    ax1.set_title('Basin Stability Correlation')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot Size Ratio
    ax2 = axes[0, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax2.plot(method_data['change_percent'], method_data['size_ratio'], 
                marker='s', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax2.set_xlabel('Change Percentage')
    ax2.set_ylabel('Size Ratio')
    ax2.set_title('Attractor Set Size Ratio')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot Functional Similarity
    ax3 = axes[1, 0]
    for method in methods:
        method_data = df[df['method'] == method]
        ax3.plot(method_data['change_percent'], method_data['functional_similarity'], 
                marker='^', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax3.set_xlabel('Change Percentage')
    ax3.set_ylabel('Functional Similarity')
    ax3.set_title('Functional Similarity Score')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot Pattern Overlap
    ax4 = axes[1, 1]
    for method in methods:
        method_data = df[df['method'] == method]
        ax4.plot(method_data['change_percent'], method_data['pattern_overlap'], 
                marker='D', linewidth=2.5, markersize=6, label=method, color=method_colors[method])
    ax4.set_xlabel('Change Percentage')
    ax4.set_ylabel('Pattern Overlap')
    ax4.set_title('Pattern Overlap Score')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'stability_functional_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create figure 4: Performance summary heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Calculate average performance for each method across all change percentages
    performance_summary = df.groupby('method')[['jaccard_similarity', 'hamming_similarity', 
                                               'composite_score', 'f1_score', 'precision', 
                                               'recall', 'functional_similarity', 'stability_correlation']].mean()
    
    # Create heatmap
    sns.heatmap(performance_summary.T, annot=True, fmt='.3f', cmap='RdYlBu_r', 
                center=0.5, ax=ax, cbar_kws={'label': 'Performance Score'})
    ax.set_title('Average Performance Summary Across All Metrics', fontsize=14, fontweight='bold')
    ax.set_xlabel('Methods')
    ax.set_ylabel('Metrics')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'performance_summary_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
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
    print(f"  - stability_functional_metrics.png")
    print(f"  - performance_summary_heatmap.png")
    if 'total_time' in df.columns:
        print(f"  - runtime_comparison.png")
    print(f"  - robustness_analysis.png")

def load_and_plot_results(csv_file_path, dataset_name='toy'):
    """
    Load results from CSV file and create all comparison plots.
    
    Parameters:
    -----------
    csv_file_path : str
        Path to the CSV file containing results
    dataset_name : str
        Name of the dataset for folder creation
    """
    # Load the data
    df = pd.read_csv(csv_file_path)
    
    # Verify required columns exist
    required_columns = ['method', 'change_percent', 'jaccard_similarity', 'hamming_similarity', 
                        'composite_score', 'precision', 'recall', 'f1_score']
    
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Warning: Missing columns: {missing_columns}")
    
    # Create all plots
    create_comparison_plots(df, dataset_name)
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print("=" * 50)
    for method in df['method'].unique():
        method_data = df[df['method'] == method]
        avg_composite = method_data['composite_score'].mean()
        avg_jaccard = method_data['jaccard_similarity'].mean()
        print(f"{method.upper():>6}: Avg Composite Score = {avg_composite:.3f}, Avg Jaccard = {avg_jaccard:.3f}")


# Example usage:
if __name__ == "__main__":
    # Replace 'your_results.csv' with the actual path to your CSV file
    csv_file_path = 'output/comparison_results_toy.csv'
    
    # Load and create plots
    load_and_plot_results(csv_file_path, dataset_name='toy')
    
    # Alternative: If you already have the DataFrame loaded
    # create_comparison_plots(your_dataframe, dataset_name='toy')