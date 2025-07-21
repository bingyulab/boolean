import pyboolnet.file_exchange as FileExchange
import pyboolnet.interaction_graphs as IG
import pyboolnet.attractors as Attractors
import pyboolnet.state_transition_graphs as STGs
from pyboolnet.repository import get_primes
from pyboolnet.file_exchange import bnet2primes
from pyboolnet.prime_implicants import create_variables
import numpy as np
import pandas as pd
from scipy.spatial.distance import jaccard, hamming
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from itertools import combinations
import warnings
from rpy2.robjects.vectors import FloatVector



def limit_float(x):
    # print(f"Limiting float: {x}, type: {type(x)}")
    if isinstance(x, FloatVector):
        x = float(x[0])
    else:
        x = float(x)
    s = str(x)
    if '.' in s and len(s.split('.')[1]) > 2:
        return round(x, 2)
    return x


class AttractorAnalysis:
    def __init__(self, ori_bnet: str, compared_bnet):
        """
        Initialize AttractorAnalysis with original and compared Boolean network files.

        Parameters:
        -----------
        ori_bnet : str
            Path to the original Boolean network file.
        compared_bnet : str or list of str
            Path(s) to the compared Boolean network file(s).
        """
        self.ori_bnet = ori_bnet
        self.compared_bnet = compared_bnet
        self.ori_bnet = ori_bnet
        self.compared_bnet = compared_bnet

    @staticmethod
    def load_bnet(bnet_file):
        return FileExchange.bnet2primes(bnet_file)

    @staticmethod
    def compute_attractors(primes):
        return Attractors.compute_attractors(primes, "synchronous")

    def get_attractors(self, bnet_file):
        primes = AttractorAnalysis.load_bnet(bnet_file)
        attrs = AttractorAnalysis.compute_attractors(primes)
        return [x['state'] for x in attrs['attractors']]
             
    def compare_attractors(self, ori_attrs, comp_attrs):
        comparator = AttractorComparison(ori_attrs, comp_attrs)
        results = comparator.comprehensive_comparison()
        return results
    
    def compare_multiple_attractors(self, ori_attrs, comp_attrs):
        results = []
        for i, reconstructed_attractors in enumerate(comp_attrs):
            comparator = AttractorComparison(ori_attrs, reconstructed_attractors)
            comparison_result = comparator.comprehensive_comparison()
            comparison_result['reconstruction_id'] = i
            results.append(comparison_result)                    
        return results
    
    def comparison(self):
        ori_attrs = self.get_attractors(self.ori_bnet)
        if isinstance(self.compared_bnet, str):
            comp_attrs = self.get_attractors(self.compared_bnet)
            results = self.compare_attractors(ori_attrs, comp_attrs)
        else:
            comp_attrs = [self.get_attractors(bnet) for bnet in self.compared_bnet]
            results = self.compare_multiple_attractors(ori_attrs, comp_attrs)
        return results
        

class AttractorComparison:
    """
    A comprehensive toolkit for comparing attractors from different network models,
    handling variable network sizes and attractor counts.
    """
    
    def __init__(self, original_attractors, reconstructed_attractors, threshold=0.9):
        """
        Initialize the comparison with original and reconstructed attractors.
        
        Parameters:
        -----------
        original_attractors : list of dict
            Original network attractors with 'dict' containing node states
        reconstructed_attractors : list of dict
            Reconstructed network attractors with 'dict' containing node states
        """
        self.original_attractors = original_attractors
        self.reconstructed_attractors = reconstructed_attractors
        self.threshold = threshold

        # Identify common nodes and project
        self.common_nodes = self._find_common_nodes()
        self.original_projected = self._project_attractors(self.original_attractors)
        self.reconstructed_projected = self._project_attractors(self.reconstructed_attractors)

        # Build binary matrices
        self.orig_vecs = self._to_matrix(self.original_projected)
        self.recon_vecs = self._to_matrix(self.reconstructed_projected)

        # Precompute similarity matrix & matching
        self._compute_matching()
        
    def _find_common_nodes(self):
        """Find nodes present in both original and reconstructed networks."""
        if not self.original_attractors or not self.reconstructed_attractors:
            return set()
        
        orig_nodes = set(self.original_attractors[0]['dict'].keys())
        recon_nodes = set(self.reconstructed_attractors[0]['dict'].keys())
        return orig_nodes.intersection(recon_nodes)
    
    def _project_attractors(self, attractors):
        """Project attractors to common node space."""
        projected = []
        for attractor in attractors:
            projected_state = {node: attractor['dict'][node] 
                             for node in self.common_nodes 
                             if node in attractor['dict']}
            projected.append(projected_state)
        return projected
    
    def _to_matrix(self, proj):
        if not proj or not self.common_nodes:
            return np.zeros((0, 0), dtype=bool)
        mat = np.array([[att[n] for n in self.common_nodes] for att in proj], dtype=bool)
        return mat
    
    def _compute_matching(self):
        n, m = len(self.orig_vecs), len(self.recon_vecs)
        size = max(n, m)
        S = np.zeros((size, size), dtype=float)

        # Fill similarity matrix
        for i in range(n):
            for j in range(m):
                # Jaccard distance on boolean vectors
                S[i, j] = 1 - jaccard(self.orig_vecs[i], self.recon_vecs[j])

        # Hungarian assignment
        cost = -S
        row, col = linear_sum_assignment(cost)
        pairs = [(i, j) for i, j in zip(row, col) if i < n and j < m]
        self.match_sims = [S[i, j] for i, j in pairs]
        
    def compute_global_jaccard(self):
        return limit_float(np.mean(self.match_sims) if self.match_sims else 0.0)

    def compute_jaccard_similarity(self):
        """
        Compute Jaccard similarity between attractor sets.
        
        Returns:
        --------
        float : Jaccard similarity coefficient (0-1)
        """
        if not self.common_nodes:
            return 0.0
        
        # Convert to binary vectors for Jaccard computation
        orig_vectors = [self._dict_to_vector(att) for att in self.original_projected]
        recon_vectors = [self._dict_to_vector(att) for att in self.reconstructed_projected]
        
        # Find best matches using minimum distance
        matches = []
        used_recon = set()
        
        for orig_vec in orig_vectors:
            best_match = None
            best_similarity = -1
            
            for i, recon_vec in enumerate(recon_vectors):
                if i in used_recon:
                    continue
                
                # Jaccard similarity = 1 - Jaccard distance
                similarity = 1 - jaccard(orig_vec, recon_vec)
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_match = i
            
            if best_match is not None:
                matches.append(best_similarity)
                used_recon.add(best_match)
            else:
                matches.append(0.0)
        
        return np.mean(matches) if matches else 0.0
    
    def compute_global_hamming(self):
        """Mean Hamming similarity over optimal one-to-one matching."""
        # Recompute with hamming if needed
        n, m = len(self.orig_vecs), len(self.recon_vecs)
        if not (n and m):
            return 0.0
        size = max(n, m)
        H = np.zeros((size, size), dtype=float)
        for i in range(n):
            for j in range(m):
                H[i, j] = 1 - hamming(self.orig_vecs[i], self.recon_vecs[j])
        cost = -H
        row, col = linear_sum_assignment(cost)
        sims = [H[i, j] for i, j in zip(row, col) if i<n and j<m]
        return limit_float(np.mean(sims) if sims else 0.0)
    
    def compute_hamming_similarity(self):
        """
        Compute normalized Hamming similarity between attractor sets.
        
        Returns:
        --------
        float : Average Hamming similarity (0-1)
        """
        if not self.common_nodes:
            return 0.0
        
        orig_vectors = [self._dict_to_vector(att) for att in self.original_projected]
        recon_vectors = [self._dict_to_vector(att) for att in self.reconstructed_projected]
        
        matches = []
        used_recon = set()
        
        for orig_vec in orig_vectors:
            best_match = None
            best_similarity = -1
            
            for i, recon_vec in enumerate(recon_vectors):
                if i in used_recon:
                    continue
                
                # Hamming similarity = 1 - normalized Hamming distance
                similarity = 1 - hamming(orig_vec, recon_vec)
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_match = i
            
            if best_match is not None:
                matches.append(best_similarity)
                used_recon.add(best_match)
            else:
                matches.append(0.0)
        
        return np.mean(matches) if matches else 0.0
    
    def compute_coverage(self):
        """
        Decoupled precision/recall with threshold.
        """
        tp = sum(1 for sim in self.match_sims if sim >= self.threshold)
        prec = tp / len(self.reconstructed_projected) if self.reconstructed_projected else 0.0
        rec  = tp / len(self.original_projected)  if self.original_projected  else 0.0
        f1   = (2*prec*rec/(prec+rec)) if (prec+rec)>0 else 0.0
        return {
            'precision': limit_float(prec),
            'recall':    limit_float(rec),
            'f1_score':  limit_float(f1),
            'true_positives': tp,
            'orig_total': len(self.original_projected),
            'recon_total': len(self.reconstructed_projected)
        }
        
    def compute_coverage_metrics(self):
        """
        Compute coverage metrics for attractor comparison.
        
        Returns:
        --------
        dict : Coverage metrics including precision, recall, and F1-score
        """
        if not self.common_nodes:
            return {'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0}
        
        orig_set = {self._dict_to_tuple(att) for att in self.original_projected}
        recon_set = {self._dict_to_tuple(att) for att in self.reconstructed_projected}
        
        intersection = orig_set.intersection(recon_set)
        
        precision = len(intersection) / len(recon_set) if recon_set else 0.0
        recall = len(intersection) / len(orig_set) if orig_set else 0.0
        f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        
        return {
            'precision': limit_float(precision),
            'recall': limit_float(recall),
            'f1_score': limit_float(f1_score),
            'exact_matches': len(intersection),
            'total_original': len(orig_set),
            'total_reconstructed': len(recon_set)
        }
    
    def compute_basin_stability_comparison(self):
        """
        Compare basin stability characteristics between attractor sets.
        
        Returns:
        --------
        dict : Basin stability comparison metrics
        """
        if not self.common_nodes:
            return {'stability_correlation': 0.0, 'size_ratio': 0.0}
        
        orig_count = len(self.original_projected)
        recon_count = len(self.reconstructed_projected)
        
        # Simple size ratio metric
        size_ratio = min(orig_count, recon_count) / max(orig_count, recon_count)
        
        # Stability correlation based on attractor diversity
        orig_diversity = self._compute_diversity(self.original_projected)
        recon_diversity = self._compute_diversity(self.reconstructed_projected)
        
        diversity_correlation = 1 - abs(orig_diversity - recon_diversity)
        
        return {
            'stability_correlation': limit_float(diversity_correlation),
            'size_ratio': limit_float(size_ratio),
            'original_count': orig_count,
            'reconstructed_count': recon_count,
            'original_diversity': limit_float(orig_diversity),
            'reconstructed_diversity': limit_float(recon_diversity)
        }
    
    def compute_functional_similarity(self):
        """
        Compute functional similarity based on active node patterns.
        
        Returns:
        --------
        dict : Functional similarity metrics
        """
        if not self.common_nodes:
            return {'functional_similarity': 0.0, 'pattern_overlap': 0.0}
        
        # Compute active node patterns
        orig_patterns = self._extract_active_patterns(self.original_projected)
        recon_patterns = self._extract_active_patterns(self.reconstructed_projected)
        
        # Pattern overlap
        common_patterns = orig_patterns.intersection(recon_patterns)
        all_patterns = orig_patterns.union(recon_patterns)
        
        pattern_overlap = len(common_patterns) / len(all_patterns) if all_patterns else 0.0
        
        # Functional similarity based on node activity correlation
        orig_activity = self._compute_node_activity(self.original_projected)
        recon_activity = self._compute_node_activity(self.reconstructed_projected)
        
        activity_correlation = np.corrcoef(
            [orig_activity.get(node, 0) for node in self.common_nodes],
            [recon_activity.get(node, 0) for node in self.common_nodes]
        )[0, 1] if len(self.common_nodes) > 1 else 0.0
        
        # Handle NaN correlation
        if np.isnan(activity_correlation):
            activity_correlation = 0.0
        
        return {
            'functional_similarity': limit_float(activity_correlation),
            'pattern_overlap': limit_float(pattern_overlap),
            'common_patterns': len(common_patterns),
            'total_patterns': len(all_patterns)
        }
    
    def comprehensive_comparison(self, Return_DF=True):
        """
        Perform comprehensive comparison using multiple metrics.
        
        Returns:
        --------
        dict : Comprehensive comparison results
        """
        results = {
            'common_nodes': len(self.common_nodes),
            # 'common_node_list': sorted(list(self.common_nodes)),
            # 'jaccard_similarity': limit_float(self.compute_jaccard_similarity()),
            # 'hamming_similarity': limit_float(self.compute_hamming_similarity()),
            # 'coverage_metrics': self.compute_coverage_metrics(),            
            'jaccard_similarity': self.compute_global_jaccard(),
            'hamming_similarity': self.compute_global_hamming(),
            'coverage_metrics': self.compute_coverage(),
            # 'basin_stability': self.compute_basin_stability_comparison(),
            # 'functional_similarity': self.compute_functional_similarity()
        }
        
        # Compute composite score
        weights = {
            'jaccard': 0.3,
            'hamming': 0.3,
            'coverage_f1': 0.4
        }
        
        composite_score = (
            weights['jaccard'] * results['jaccard_similarity'] +
            weights['hamming'] * results['hamming_similarity'] +
            weights['coverage_f1'] * results['coverage_metrics']['f1_score'] 
        )
        
        results['composite_score'] = limit_float(composite_score)
        
        return self.to_dataframe(results) if Return_DF else results
    
    def to_dataframe(self, results):
        flat = results.copy()
        if 'coverage_metrics' in flat:
            flat['coverage_metrics'] = flat.pop('coverage_metrics')
        if 'basin_stability' in flat:
            flat['basin_stability'] = flat.pop('basin_stability')
        if 'functional_similarity' in flat:
            flat['functional_similarity'] = flat.pop('functional_similarity')

        # Create DataFrame
        df = pd.DataFrame([flat])
        return df
    
    def _dict_to_vector(self, attractor_dict):
        """Convert attractor dictionary to binary vector."""
        return np.array([attractor_dict.get(node, 0) for node in sorted(self.common_nodes)])
    
    def _dict_to_tuple(self, attractor_dict):
        """Convert attractor dictionary to tuple for set operations."""
        return tuple(attractor_dict.get(node, 0) for node in sorted(self.common_nodes))
    
    def _compute_diversity(self, attractors):
        """Compute diversity measure for attractor set."""
        if len(attractors) <= 1:
            return 0.0
        
        vectors = [self._dict_to_vector(att) for att in attractors]
        total_distance = 0
        count = 0
        
        for i, j in combinations(range(len(vectors)), 2):
            total_distance += hamming(vectors[i], vectors[j])
            count += 1
        
        return total_distance / count if count > 0 else 0.0
    
    def _extract_active_patterns(self, attractors):
        """Extract active node patterns from attractors."""
        patterns = set()
        for attractor in attractors:
            active_nodes = frozenset(node for node, state in attractor.items() if state == 1)
            patterns.add(active_nodes)
        return patterns
    
    def _compute_node_activity(self, attractors):
        """Compute average activity level for each node."""
        activity = {}
        for node in self.common_nodes:
            activity[node] = np.mean([att.get(node, 0) for att in attractors])
        return activity


# Example usage and demonstration
def demonstrate_attractor_comparison():
    """Demonstrate the attractor comparison functionality."""
    
    # Example original attractors (from your data)
    original_attractors = [
        {'dict': {'Akt': 0, 'EGF': 0, 'Erk': 1, 'Hsp27': 1, 'Jnk': 0, 'Mek': 1, 
                 'NFkB': 0, 'PI3K': 0, 'Raf': 0, 'TNFa': 0, 'cJun': 0, 'p90RSK': 1}},
        {'dict': {'Akt': 1, 'EGF': 0, 'Erk': 0, 'Hsp27': 1, 'Jnk': 1, 'Mek': 0, 
                 'NFkB': 1, 'PI3K': 1, 'Raf': 0, 'TNFa': 1, 'cJun': 1, 'p90RSK': 0}},
        {'dict': {'Akt': 1, 'EGF': 1, 'Erk': 1, 'Hsp27': 1, 'Jnk': 0, 'Mek': 1, 
                 'NFkB': 0, 'PI3K': 1, 'Raf': 1, 'TNFa': 0, 'cJun': 0, 'p90RSK': 1}},
        {'dict': {'Akt': 1, 'EGF': 1, 'Erk': 1, 'Hsp27': 1, 'Jnk': 1, 'Mek': 1, 
                 'NFkB': 1, 'PI3K': 1, 'Raf': 1, 'TNFa': 1, 'cJun': 1, 'p90RSK': 1}}
    ]
    
    # Example reconstructed attractors (smaller network, different count)
    reconstructed_attractors = [
        {'dict': {'Akt': 1, 'EGF': 1, 'Erk': 1, 'Hsp27': 1, 'Mek': 1, 'p90RSK': 1}},
        {'dict': {'Akt': 0, 'EGF': 0, 'Erk': 0, 'Hsp27': 0, 'Mek': 0, 'p90RSK': 0}},
        {'dict': {'Akt': 1, 'EGF': 0, 'Erk': 1, 'Hsp27': 1, 'Mek': 1, 'p90RSK': 1}}
    ]
    
    # Perform comparison
    comparator = AttractorComparison(original_attractors, reconstructed_attractors)
    results = comparator.comprehensive_comparison()
    
    # Display results
    print("Attractor Comparison Results:")
    print(f"Common nodes: {results['common_nodes']}")
    # print(f"Common node list: {results['common_node_list']}")
    print(f"Jaccard similarity: {results['jaccard_similarity'].iloc[0]:.3f}")
    print(f"Hamming similarity: {results['hamming_similarity'].iloc[0]:.3f}")
    print(f"Coverage F1-score: {results['coverage_metrics'].iloc[0]['f1_score']:.3f}")
    # print(f"Functional similarity: {results['functional_similarity']['functional_similarity']:.3f}")
    print(f"Composite score: {results['composite_score'].iloc[0]:.3f}")
    
    return results

# Utility functions for batch processing
def compare_multiple_reconstructions(original_attractors, reconstruction_list):
    """
    Compare original attractors with multiple reconstructions.
    
    Parameters:
    -----------
    original_attractors : list of dict
        Original network attractors
    reconstruction_list : list of list of dict
        List of reconstructed attractor sets
    
    Returns:
    --------
    list : Comparison results for each reconstruction
    """
    results = []
    for i, reconstructed_attractors in enumerate(reconstruction_list):
        comparator = AttractorComparison(original_attractors, reconstructed_attractors)
        comparison_result = comparator.comprehensive_comparison()
        comparison_result['reconstruction_id'] = i
        results.append(comparison_result)
    
    return results

def create_comparison_summary(comparison_results):
    """
    Create summary statistics from multiple comparisons.
    
    Parameters:
    -----------
    comparison_results : list of dict
        Results from multiple attractor comparisons
    
    Returns:
    --------
    dict : Summary statistics
    """
    if not comparison_results:
        return {}
    
    metrics = ['jaccard_similarity', 'hamming_similarity', 'composite_score']
    summary = {}
    
    for metric in metrics:
        values = [result[metric] for result in comparison_results]
        summary[metric] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'median': np.median(values)
        }
    
    # Coverage metrics summary
    f1_scores = [result['coverage_metrics']['f1_score'] for result in comparison_results]
    summary['coverage_f1'] = {
        'mean': np.mean(f1_scores),
        'std': np.std(f1_scores),
        'min': np.min(f1_scores),
        'max': np.max(f1_scores),
        'median': np.median(f1_scores)
    }
    
    return summary

if __name__ == "__main__":
    # Run demonstration
    demonstrate_attractor_comparison()
