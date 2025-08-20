import pyboolnet.file_exchange as FileExchange
import pyboolnet.attractors as Attractors
import pyboolnet.basins_of_attraction as Basins
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from rpy2.robjects.vectors import FloatVector
from scipy.stats import entropy


def limit_float(x, nbit=2):
    # print(f"Limiting float: {x}, type: {type(x)}")
    if isinstance(x, FloatVector):
        x = float(x[0])
    else:
        x = float(x)
    s = str(x)
    if '.' in s and len(s.split('.')[1]) > nbit:
        return round(x, nbit)
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
        import logging 
        logging.disable(logging.INFO)
        result = Attractors.compute_attractors(primes, "synchronous") 
        logging.disable(logging.NOTSET)
        return result

    def get_attractors(self, bnet_file):
        primes = AttractorAnalysis.load_bnet(bnet_file)
        attrs = AttractorAnalysis.compute_attractors(primes)
        return primes, [x['state'] for x in attrs['attractors']]
    
    @staticmethod
    def get_basin_sizes(primes, state):         
        import logging 
        logging.disable(logging.INFO)
        weak = Basins.weak_basin(primes, "asynchronous", state)
        logging.disable(logging.NOTSET)
        return weak.get('perc', 0.0)
    
    def compare_attractors(self, ori_primes, ori_attrs, com_primes, comp_attrs):
        comparator = AttractorComparison(ori_primes, ori_attrs, com_primes, comp_attrs)
        results = comparator.comprehensive_comparison()
        return results
    
    def compare_multiple_attractors(self, ori_primes, ori_attrs, com_primes, comp_attrs):
        results = []
        for i, reconstructed_attractors in enumerate(comp_attrs):
            comparator = AttractorComparison(ori_primes, ori_attrs, com_primes[i], reconstructed_attractors)
            comparison_result = comparator.comprehensive_comparison()
            comparison_result['reconstruction_id'] = i
            results.append(comparison_result)                    
        return results
    
    def comparison(self):
        ori_primes, ori_attrs = self.get_attractors(self.ori_bnet)
        if isinstance(self.compared_bnet, str):
            comp_primes, comp_attrs = self.get_attractors(self.compared_bnet)
            results = self.compare_attractors(ori_primes, ori_attrs, comp_primes, comp_attrs)
        else:
            tmp_res = [self.get_attractors(bnet) for bnet in self.compared_bnet]
            com_primes = [x[0] for x in tmp_res]
            comp_attrs = [x[1] for x in tmp_res]
            results = self.compare_multiple_attractors(ori_primes, ori_attrs, com_primes, comp_attrs)
        return results
            

class AttractorComparison:
    """
    A comprehensive toolkit for comparing attractors from different network models,
    handling variable network sizes and attractor counts.
    """

    def __init__(self, ori_primes, original_attractors, recon_primes, 
                 reconstructed_attractors, basin_weight=0.3):
        self.ori_primes = ori_primes
        self.recon_primes = recon_primes
        self.basin_weight = basin_weight  # Weight for basin size consideration        
        
        self.recon_total = len(reconstructed_attractors)
        self.original_attractors = original_attractors
        self.reconstructed_attractors = self._select_attractors_for_comparison(
            reconstructed_attractors, len(original_attractors)
        )
        
        # Get node sets
        self.orig_nodes = self._get_nodes_from_attractors(self.original_attractors)
        self.recon_nodes = self._get_nodes_from_attractors(self.reconstructed_attractors)
        self.common_nodes = self.orig_nodes.intersection(self.recon_nodes)
        self.all_nodes = sorted(list(self.orig_nodes.union(self.recon_nodes)))
        
        # Build enhanced representations
        self._build_enhanced_representations()
        
        # Compute basin sizes for weighting
        self.orig_basins = self._compute_basin_sizes(ori_primes, self.original_attractors)
        self.recon_basins = self._compute_basin_sizes(recon_primes, self.reconstructed_attractors)
            
    def _compute_basin_sizes(self, primes, attractors):
        weight = [AttractorAnalysis.get_basin_sizes(primes, att['str']) for att in attractors]        
        weights = np.array(weight)
        # Normalize to prevent scale issues
        return weights / np.sum(weights) if np.sum(weights) > 0 else weights
    
    def _select_attractors_for_comparison(self, attractors, max_count):
        if len(attractors) <= max_count:
            return attractors
        # Select top K attractors by basin size
        basin_sizes = self._compute_basin_sizes(self.recon_primes, attractors)
        #  Normalize basin sizes to [0,1] to prevent overflow
        basin_sizes = basin_sizes / np.max(basin_sizes) if np.max(basin_sizes) > 0 else basin_sizes
        # Get indices of top K attractors
        top_indices = np.argsort(basin_sizes)[::-1][:max_count]
        return [attractors[i] for i in top_indices]            
    
    def _get_nodes_from_attractors(self, attractors):
        """Get all unique nodes from attractors."""
        if not attractors:
            return set()
        nodes = set()
        for att in attractors:
            nodes.update(att['dict'].keys())
        return nodes
    
    def _build_enhanced_representations(self):
        """
        Build enhanced representations that handle missing nodes properly.
        Uses three-state representation: 0, 1, NaN (missing)
        """
        self.orig_enhanced = []
        self.recon_enhanced = []
        
        # Build for original attractors
        for att in self.original_attractors:
            vector = []
            for node in self.all_nodes:
                if node in att['dict']:
                    vector.append(float(att['dict'][node]))
                else:
                    vector.append(np.nan)  # Missing node
            self.orig_enhanced.append(vector)
        
        # Build for reconstructed attractors
        for att in self.reconstructed_attractors:
            vector = []
            for node in self.all_nodes:
                if node in att['dict']:
                    vector.append(float(att['dict'][node]))
                else:
                    vector.append(np.nan)  # Missing node
            self.recon_enhanced.append(vector)
        
        self.orig_enhanced = np.array(self.orig_enhanced)
        self.recon_enhanced = np.array(self.recon_enhanced)
    
    def _compute_pairwise_similarity(self, vec1, vec2, similarity_type='jaccard'):
        """
        Compute similarity between two vectors handling NaN values properly.
        """
        # Find positions where both vectors have defined values
        mask = ~(np.isnan(vec1) | np.isnan(vec2))
        
        if not np.any(mask):
            return 0.0  # No common nodes
        
        v1_masked = vec1[mask]
        v2_masked = vec2[mask]
        
        if similarity_type == 'jaccard':
            # Modified Jaccard for missing nodes
            intersection = np.sum(v1_masked == v2_masked)
            union = len(v1_masked)  # All compared positions
            return intersection / union if union > 0 else 0.0
            
        elif similarity_type == 'hamming':
            # Hamming similarity
            matches = np.sum(v1_masked == v2_masked)
            total = len(v1_masked)
            return matches / total if total > 0 else 0.0

    def _compute_similarity_matrix(self, similarity_type='jaccard'):
        """Compute similarity matrix between all pairs of attractors."""
        n_orig = len(self.orig_enhanced)
        n_recon = len(self.recon_enhanced)
        
        if n_orig == 0 or n_recon == 0:
            return np.array([]).reshape(n_orig, n_recon)
        
        sim_matrix = np.zeros((n_orig, n_recon))
        
        for i in range(n_orig):
            for j in range(n_recon):
                sim_matrix[i, j] = self._compute_pairwise_similarity(
                    self.orig_enhanced[i], self.recon_enhanced[j], similarity_type)
        
        return sim_matrix

    def compute_optimal_matching_metrics(self):
        """
        Compute metrics based on optimal bipartite matching.
        This addresses the third problem by finding best matches first.
        """
        # Compute similarity matrices for different metrics
        jaccard_matrix = self._compute_similarity_matrix('jaccard')
        hamming_matrix = self._compute_similarity_matrix('hamming')
        
        n_orig = len(self.original_attractors)
        n_recon = len(self.reconstructed_attractors)
        
        if n_orig == 0 or n_recon == 0:
            return self._empty_results()
        
        # Use Jaccard for matching (most robust for categorical data)
        cost_matrix = 1 - jaccard_matrix
        
        # Handle rectangular matrices for Hungarian algorithm
        if n_orig != n_recon:
            # Pad with high costs to handle different sizes
            max_dim = max(n_orig, n_recon)
            padded_cost = np.full((max_dim, max_dim), 2.0)  # Cost > 1
            padded_cost[:n_orig, :n_recon] = cost_matrix
            row_ind, col_ind = linear_sum_assignment(padded_cost)
            
            # Filter out padding matches
            valid_matches = [(i, j) for i, j in zip(row_ind, col_ind) 
                           if i < n_orig and j < n_recon]
        else:
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            valid_matches = list(zip(row_ind, col_ind))
        
        # Compute metrics for matched pairs
        matched_jaccard = []
        matched_hamming = []
        
        for i, j in valid_matches:
            matched_jaccard.append(jaccard_matrix[i, j])
            matched_hamming.append(hamming_matrix[i, j])
        
        true_positives = 0
        false_positives = 0
        false_negatives = 0

        for i, j in valid_matches:
            vec_true = self.orig_enhanced[i]
            vec_pred = self.recon_enhanced[j]
            
            # Iterate over all nodes in the union space
            for node_idx in range(len(self.all_nodes)):
                true_val = vec_true[node_idx]
                pred_val = vec_pred[node_idx]
                
                # Valid comparison is possible only if the original value is not NaN
                if not np.isnan(true_val):
                    # True Positive: Both are 1 and match
                    if true_val == 1 and pred_val == 1:
                        true_positives += 1
                    # False Positive: Original is 0, but prediction is 1
                    elif true_val == 0 and pred_val == 1:
                        false_positives += 1
                    # False Negative: Original is 1, but prediction is 0 or NaN (missing)
                    elif true_val == 1 and (pred_val == 0 or np.isnan(pred_val)):
                        false_negatives += 1

        # Calculate precision, recall, and F1 score from the counts
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0.0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0.0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        return {
            'jaccard': np.mean(matched_jaccard) if matched_jaccard else 0.0,
            'hamming': np.mean(matched_hamming) if matched_hamming else 0.0,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'total_matches': len(valid_matches),
            'orig_count': n_orig,
            'recon_count': n_recon
        }

    def compute_structural_similarity(self):
        """
        Compute structural similarity accounting for attractor count differences.
        This addresses the second problem.
        """
        n_orig = len(self.original_attractors)
        n_recon = len(self.reconstructed_attractors)
        
        # Count penalty for different numbers of attractors
        count_penalty = abs(n_orig - n_recon) / max(n_orig, n_recon, 1)
        
        # Basin size distribution similarity (if available)
        basin_similarity = 0.0
        if len(self.orig_basins) > 0 and len(self.recon_basins) > 0:
            # Compare basin size distributions using KL divergence
            orig_dist = self.orig_basins + 1e-10  # Add small value to avoid log(0)
            recon_dist = self.recon_basins + 1e-10
            
            # Normalize to probabilities
            orig_dist = orig_dist / np.sum(orig_dist)
            recon_dist = recon_dist / np.sum(recon_dist)
            
            # Pad shorter distribution
            max_len = max(len(orig_dist), len(recon_dist))
            if len(orig_dist) < max_len:
                orig_dist = np.pad(orig_dist, (0, max_len - len(orig_dist)), 'constant', constant_values=1e-10)
            if len(recon_dist) < max_len:
                recon_dist = np.pad(recon_dist, (0, max_len - len(recon_dist)), 'constant', constant_values=1e-10)
            
            # Jensen-Shannon divergence (symmetric and bounded)
            m = 0.5 * (orig_dist + recon_dist)
            js_div = 0.5 * (entropy(orig_dist, m) + entropy(recon_dist, m))
            basin_similarity = 1 - js_div  # Convert divergence to similarity
        
        structural_score = (1 - count_penalty) * (1 - self.basin_weight) + basin_similarity * self.basin_weight
        
        return {
            'structural_score': max(0, structural_score),  # Ensure non-negative
            'count_penalty': count_penalty,
            'basin_similarity': basin_similarity,
            'orig_attractor_count': n_orig,
            'recon_attractor_count': n_recon,            
            'recon_total': self.recon_total
        }
    
    def _empty_results(self):
        """Return empty results when comparison is not possible."""
        return {
            'jaccard': 0.0, 'hamming': 0.0, 
            'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0,
            'high_quality_matches': 0, 'total_matches': 0,
            'orig_count': len(self.original_attractors),
            'recon_count': len(self.reconstructed_attractors),
            'recon_total': self.recon_total
        }
    
    def comprehensive_comparison(self, return_df=True):
        """
        Perform comprehensive comparison with improved metrics.
        """
        # Optimal matching metrics
        matching_results = self.compute_optimal_matching_metrics()
        
        # Structural similarity
        structural_results = self.compute_structural_similarity()
        
        # Node coverage analysis
        node_results = {
            'common_nodes': len(self.common_nodes),
            'orig_nodes': len(self.orig_nodes),
            'recon_nodes': len(self.recon_nodes),
            'node_coverage': len(self.common_nodes) / len(self.orig_nodes) if self.orig_nodes else 0.0
        }
        
        # Combine all results
        results = {**matching_results, **structural_results, **node_results}
        
        # Compute composite score with balanced weighting
        composite_score = (
            0.3 * results['jaccard'] +
            0.2 * results['hamming'] +
            0.3 * results['f1_score'] +
            0.2 * results['structural_score']
        )
        
        results['composite_score'] = composite_score
        
        if return_df:
            return pd.DataFrame([results])
        return results     
    
    
if __name__ == "__main__":
    ori_primes = "data/ToyModel/ToyModel.bnet"
    # compared_bnet = "data/ToyModel/ToyModel.bnet"
    # compared_bnet = "output/cellnopt/ToyModel/80_Modified/ga/OPT_ToyModel.bnet"
    compared_bnet = "output/cellnopt/ToyModel/90_Modified/ga/OPT_ToyModel.bnet"
    analysis = AttractorAnalysis(ori_primes, compared_bnet)
    results = analysis.comparison()
    print(results)