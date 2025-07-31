import pyboolnet.file_exchange as FileExchange
import pyboolnet.interaction_graphs as IG
import pyboolnet.attractors as Attractors
import pyboolnet.basins_of_attraction as Basins
import pyboolnet.state_transition_graphs as STGs
import numpy as np
import pandas as pd
from scipy.spatial.distance import jaccard, hamming
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from itertools import combinations
import warnings
from rpy2.robjects.vectors import FloatVector


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

    def __init__(self, ori_primes, original_attractors, recon_primes, reconstructed_attractors, threshold=0.9):
        self.ori_primes = ori_primes
        self.com_primes = recon_primes
        
        self.original_attractors = original_attractors
        self.reconstructed_attractors = reconstructed_attractors
        self.threshold = threshold

        # Identify common nodes and project                
        self.common_nodes = self._find_common_nodes()
        # Build union of all nodes from both networks
        ori_nodes = set(self.original_attractors[0]['dict'].keys())
        recon_nodes = set(self.reconstructed_attractors[0]['dict'].keys())
        self.all_nodes = sorted(list(ori_nodes.union(recon_nodes)))
    
        self.orig_vecs, self.pres_ori = self._build_matrices(self.original_attractors)
        self.recon_vecs, self.pres_recon = self._build_matrices(self.reconstructed_attractors)

        self.ori_basin_sizes = self._compute_basin_sizes(self.ori_primes, self.original_attractors)
        self.recon_basin_sizes = self._compute_basin_sizes(self.com_primes, self.reconstructed_attractors)
        
        b_ori_norm = self.ori_basin_sizes / (self.ori_basin_sizes.sum() or 1)
        b_recon_norm = self.recon_basin_sizes / (self.recon_basin_sizes.sum() or 1)
        # Geometric mean weights
        self.W = np.sqrt(b_ori_norm.reshape(-1, 1) * b_recon_norm.reshape(1, -1))                           # shape (n,m)

        # Precompute enhanced similarity matrix & matching
        self._compute_matching()
        
    def _find_common_nodes(self):
        """Find nodes present in both original and reconstructed networks."""
        if not self.original_attractors or not self.reconstructed_attractors:
            return set()
        
        orig_nodes = set(self.original_attractors[0]['dict'].keys())
        recon_nodes = set(self.reconstructed_attractors[0]['dict'].keys())
        return orig_nodes.intersection(recon_nodes)    
    
    def _build_matrices(self, attractors):
        """
        For each attractor, build:
          - vec: boolean vector over all_nodes (True if node=1)
          - pres: boolean mask (True if node is defined in the attractor)
        """
        vecs = []
        pres = []
        for att in attractors:
            d = att['dict']
            vecs.append([bool(d.get(node, False)) for node in self.all_nodes])
            pres.append([node in d for node in self.all_nodes])
        return np.array(vecs, dtype=bool), np.array(pres, dtype=bool)
    
    def _compute_basin_sizes(self, primes, attractors):
        weight = [AttractorAnalysis.get_basin_sizes(primes, att['str']) for att in attractors]
        return np.array(weight)

    def _compute_matching(self):
        n, m = len(self.orig_vecs), len(self.recon_vecs)
        if n == 0 or m == 0:
            self.match_sims = []
            return
        
        # Compute masked Jaccard for each pair
        sim_matrix = np.zeros((n, m), dtype=float)
        for i in range(n):
            for j in range(m):
                mask = self.pres_ori[i] & self.pres_recon[j]
                if not mask.any():
                    sim_matrix[i, j] = 0.0
                else:
                    a = self.orig_vecs[i]
                    b = self.recon_vecs[j]
                    inter = np.sum((a & b) & mask)
                    uni   = np.sum((a | b) & mask)
                    sim_matrix[i, j] = inter / uni if uni > 0 else 0.0

        # Hungarian assignment
        cost = -(self.W * sim_matrix)
        row, col = linear_sum_assignment(cost)
        pairs = [(i, j) for i, j in zip(row, col) if i < n and j < m]
        self.match_sims = [sim_matrix[i, j] for i, j in pairs]

    def compute_global_jaccard(self):
        # larger values indeed mean better matches
        return limit_float(np.mean(self.match_sims)  if self.match_sims else 0.0)

    def compute_global_hamming(self):
        """Mean Hamming similarity over optimal one-to-one matching."""
        # Larger values indicate better matches
        # Recompute with hamming if needed
        n, m = len(self.orig_vecs), len(self.recon_vecs)
        if len(self.orig_vecs) == 0 or len(self.recon_vecs) == 0:
            return 0.0        
        
        # Compute masked Hamming for each pair
        # D_{ij} = sum_k |a_i[k] - b_j[k]| over only coords where both define the node.    
        dist_matrix = np.zeros((n, m), dtype=float)
        for i in range(n):
            for j in range(m):
                mask = self.pres_ori[i] & self.pres_recon[j]
                if not mask.any():
                    dist_matrix[i, j] = 0.0
                else:
                    a = self.orig_vecs[i, mask].astype(int)
                    b = self.recon_vecs[j, mask].astype(int)
                    # Hamming distance = sum of mismatches
                    dist_matrix[i, j] = np.sum(np.abs(a - b))
                
        # For Hungarian algorithm, we use distances directly as costs
        # (no need to negate since we want minimum distance)
        row, col = linear_sum_assignment(self.W * dist_matrix)
        
        # Extract the matched similarities
        match_sims = [dist_matrix[i, j] for i, j in zip(row, col)
                      if i < len(self.orig_vecs) and j < len(self.recon_vecs)]
        
        if not match_sims:
            return 0.0

        # Convert mean Hamming distance back into similarity [0,1]:
        # max distance per pair = #masked_coords, but these vary per pair.
        # We approximate by dividing by the average masked-length:
        lengths = []
        for i, j in zip(row, col):
            mask = self.pres_ori[i] & self.pres_recon[j]
            lengths.append(np.sum(mask))
        avg_length = np.mean(lengths) if lengths else 1.0

        mean_dist = np.mean(match_sims)
        similarity = 1.0 - (mean_dist / avg_length)
        return limit_float(similarity)

    def compute_coverage(self):
        """
        Decoupled precision/recall with threshold.
        """
        tp = sum(1 for sim in self.match_sims if sim >= self.threshold)
        prec = tp / len(self.reconstructed_attractors) if self.reconstructed_attractors else 0.0
        rec  = tp / len(self.original_attractors)  if self.original_attractors  else 0.0
        f1   = (2*prec*rec/(prec+rec)) if (prec+rec)>0 else 0.0
        return {
            'precision': limit_float(prec),
            'recall':    limit_float(rec),
            'f1_score':  limit_float(f1),
            'true_positives': tp,
            'orig_total': len(self.original_attractors),
            'recon_total': len(self.reconstructed_attractors)
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
            'jaccard_similarity': self.compute_global_jaccard(),
            'hamming_similarity': self.compute_global_hamming(),
            'coverage_metrics': self.compute_coverage(),     
        }
        
        # Compute composite score
        weights = {
            'jaccard': 0.4,
            'hamming': 0.4,
            'coverage_f1': 0.2
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
            cov = flat.pop('coverage_metrics')
            for k, v in cov.items():
                flat[k] = v
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

    
    
if __name__ == "__main__":
    ori_primes = "data/ToyModel/ToyModel.bnet"
    compared_bnet = "data/ToyModel/ToyModel.bnet"
    analysis = AttractorAnalysis(ori_primes, compared_bnet)
    results = analysis.comparison()
    print(results)