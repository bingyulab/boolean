import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats


class HybridBayesianPKN:
    def __init__(self, alpha=0.05, biological_prior_weight=1.0, structure_penalty=0.1):
        """
        Hybrid Bayesian Network Learning for PKN Discovery
        
        Parameters:
        - alpha: significance level for statistical tests
        - biological_prior_weight: weight for biological prior knowledge
        - structure_penalty: penalty for complex structures (L1-like regularization)
        """
        self.alpha = alpha
        self.biological_prior_weight = biological_prior_weight
        self.structure_penalty = structure_penalty
        self.learned_network = {}
        self.edge_scores = {}
        
    def parse_midas_data(self, medias_data):
        """Parse MIDAS format data"""
        print(f"Parsing MIDAS data...{medias_data}")
        df = pd.read_csv(medias_data)
        df = df.apply(pd.to_numeric, errors='coerce')
        
        # Separate treatments, drug applications, and measurements
        tr_cols = [col for col in df.columns if col.startswith('TR:') and col != 'TR:mock:CellLine']
        da_cols = [col for col in df.columns if col.startswith('DA:')]
        dv_cols = [col for col in df.columns if col.startswith('DV:')]
        
        return df, tr_cols, da_cols, dv_cols
    
    def parse_pkn_prior(self, pkn_data):
        """Parse PKN prior knowledge in CSV format"""
        with open(pkn_data, 'r') as f:
            lines = f.read().strip().split('\n')
            prior_network = {}
            
            for line in lines[1:]:  # Skip header
                parts = line.split(' , ')
                if len(parts) == 2:
                    target, factors = parts[0].strip(), parts[1].strip()
                    # Parse boolean logic (simplified)
                    if '|' in factors:
                        # OR logic - multiple possible parents
                        parents = [f.strip().replace('!', '') for f in factors.split('|')]
                    else:
                        parents = [factors.replace('!', '').strip()]
                    prior_network[target] = parents
                
        return prior_network
    
    def normalize_node_name(self, node_name):
        """Normalize node names - remove prefixes and inhibitor 'i' suffix"""
        # Remove column prefixes
        clean_name = node_name.replace('TR:', '').replace('DA:', '').replace('DV:', '')
        # Remove inhibitor suffix 'i' (e.g., 'Rafi' -> 'Raf', 'PI3Ki' -> 'PI3K')
        if clean_name.endswith('i') and len(clean_name) > 1:
            clean_name = clean_name[:-1]
        return clean_name
    
    def calculate_biological_prior(self, source, target, prior_network):
        """
        Biological Prior Layer
        Calculate biological prior score for an edge
        """
        norm_source = self.normalize_node_name(source)
        norm_target = self.normalize_node_name(target)        

        if norm_target in prior_network:
            if norm_source in prior_network[norm_target]:
                return 2.0  # Strong prior support
            else:
                # Check if source is in same pathway (transitivity)
                for parent in prior_network[norm_target]:
                    if parent in prior_network and norm_source in prior_network.get(parent, []):
                        return 1.0  # Moderate prior support (indirect)# Check reverse direction for bidirectional relationships
        if norm_source in prior_network and norm_target in prior_network[norm_source]:
            return 1.5  # Moderate-strong support (reverse direction)            
        return 0.1  # Weak prior support
    
    def test_conditional_independence(self, data, x, y, z_set=None):
        """
        Statistical Evidence Layer: 
        Test conditional independence using partial correlation
        """
        if z_set is None or len(z_set) == 0:
            # Simple correlation test
            corr, p_value = stats.pearsonr(data[x], data[y])
            return p_value > self.alpha, abs(corr), p_value
        
        # Partial correlation test
        try:
            # Prepare data
            X_data = data[[x] + list(z_set)].values
            y_data = data[y].values
            
            # Fit regression models
            reg_xy = LinearRegression().fit(data[list(z_set)].values, data[x].values)
            reg_zy = LinearRegression().fit(data[list(z_set)].values, y_data)
            
            # Get residuals
            x_resid = data[x].values - reg_xy.predict(data[list(z_set)].values)
            y_resid = y_data - reg_zy.predict(data[list(z_set)].values)
            
            # Test correlation of residuals
            if np.std(x_resid) > 1e-10 and np.std(y_resid) > 1e-10:
                corr, p_value = stats.pearsonr(x_resid, y_resid)
                return p_value > self.alpha, abs(corr), p_value
            else:
                return True, 0.0, 1.0
                
        except:
            return True, 0.0, 1.0
    
    def calculate_treatment_effect(self, data, treatment, outcome, confounders=None):
        """
        Causal Evidence Layer
        Calculate average treatment effect
        """
        try:
            if confounders is None:
                confounders = []
            
            # Prepare features
            features = confounders if confounders else ['intercept']
            if features == ['intercept']:
                X = np.ones((len(data), 1))
            else:
                X = data[confounders].values
                
            # Split by treatment
            treated = data[data[treatment] == 1]
            control = data[data[treatment] == 0]
            
            if len(treated) == 0 or len(control) == 0:
                return 0.0, 1.0
            
            # Calculate simple difference for now
            treated_outcome = treated[outcome].mean()
            control_outcome = control[outcome].mean()
            ate = treated_outcome - control_outcome
            
            # Simple t-test for significance
            _, p_value = stats.ttest_ind(treated[outcome], control[outcome])
            
            return ate, p_value
            
        except:
            return 0.0, 1.0
    
    def score_edge(self, data, source, target, prior_network, all_nodes):
        """Score a potential edge using hybrid approach"""
        
        # 1. Statistical evidence
        is_independent, corr_strength, p_value = self.test_conditional_independence(
            data, source, target
        )
        statistical_score = corr_strength if not is_independent else 0.0
        
        # 2. Causal evidence (treatment effect if applicable)
        causal_score = 0.0
        if source.startswith('TR:'):
            ate, ate_p = self.calculate_treatment_effect(data, source, target)
            if ate_p < self.alpha:
                causal_score = abs(ate)
        
        # 3. Biological prior
        biological_score = self.calculate_biological_prior(source, target, prior_network)
        
        # 4. Combined score
        combined_score = (
            statistical_score + 
            causal_score + 
            self.biological_prior_weight * biological_score - 
            self.structure_penalty
        )
        
        return combined_score, {
            'statistical': statistical_score,
            'causal': causal_score, 
            'biological': biological_score,
            'p_value': p_value
        }
    
    def learn_network_structure(self, data, tr_cols, da_cols, dv_cols, prior_network):
        """Learn network structure using hybrid approach"""
        
        # All potential nodes
        all_nodes = tr_cols + da_cols + dv_cols
        
        # Only learn edges to measurement nodes (DA, DV)
        target_nodes = da_cols + dv_cols
        source_nodes = tr_cols + da_cols  # Treatments and intermediate measurements
        
        edge_candidates = []
        
        # Score all potential edges
        for source in source_nodes:
            for target in target_nodes:
                if source != target:
                    score, details = self.score_edge(data, source, target, prior_network, all_nodes)
                    
                    edge_candidates.append({
                        'source': source,
                        'target': target,
                        'score': score,
                        'details': details
                    })
        
        # Sort by score and select edges
        edge_candidates.sort(key=lambda x: x['score'], reverse=True)
        
        # Select edges above threshold or top K
        selected_edges = []
        # Strategy 1: Adaptive threshold based on score distribution
        scores = [edge['score'] for edge in edge_candidates]
        if len(scores) > 0:
            score_mean = np.mean(scores)
            score_std = np.std(scores) if len(scores) > 1 else 0
            adaptive_threshold = max(0.5, score_mean + 0.5 * score_std)
        else:
            adaptive_threshold = 0.5
        
        # Strategy 2: Limit edges per target (prevent hub nodes)
        edges_per_target = {}
        max_parents_per_node = 3  # Biological constraint
        
        for edge in edge_candidates:
            target = edge['target']
            
            # Multiple selection criteria
            criteria_met = (
                edge['score'] > adaptive_threshold and
                edge['details']['p_value'] < self.alpha and
                (edge['details']['statistical'] > 0.3 or edge['details']['causal'] > 0.1 or edge['details']['biological'] > 1.0) and
                edges_per_target.get(target, 0) < max_parents_per_node
            )
            
            if criteria_met:
                selected_edges.append(edge)
                edges_per_target[target] = edges_per_target.get(target, 0) + 1
        
        # Strategy 3: If still too many edges, keep only top K per target
        if len(selected_edges) > len(target_nodes) * 2:  # Max 2 edges per target on average
            final_edges = []
            target_edge_counts = {}
            
            for edge in sorted(selected_edges, key=lambda x: x['score'], reverse=True):
                target = edge['target']
                if target_edge_counts.get(target, 0) < 2:  # Max 2 parents per node
                    final_edges.append(edge)
                    target_edge_counts[target] = target_edge_counts.get(target, 0) + 1
                    
            selected_edges = final_edges
                
        print(f"Selected {len(selected_edges)} edges from {len(edge_candidates)} candidates")
        print(f"Adaptive threshold: {adaptive_threshold:.3f}")
        
        return selected_edges
    
    def infer_edge_sign(self, data, source, target):
        """Infer edge sign (positive/negative regulation)"""
        try:
            # Calculate correlation
            corr, _ = stats.pearsonr(data[source], data[target])
            return 1 if corr > 0 else -1
        except:
            return 1
    
    def convert_to_bnet(self, output_file=None):
        """Convert learned network to BNET format"""
        if not hasattr(self, 'learned_edges'):
            raise ValueError("Must fit model first")            
        # Collect all nodes
        all_nodes = set()
        node_parents = {}        
        for edge in self.learned_edges:
            source = self.normalize_node_name(edge['source'])
            target = self.normalize_node_name(edge['target'])            
            all_nodes.add(source)
            all_nodes.add(target)            
            if target not in node_parents:
                node_parents[target] = []            
            # Store edge with sign information
            sign = self.infer_edge_sign(self.data, edge['source'], edge['target'])
            node_parents[target].append((source, sign))        
        # Generate BNET format
        bnet_lines = []        
        for node in sorted(all_nodes):
            if node in node_parents:
                parents = node_parents[node]                
                if len(parents) == 1:
                    # Single parent
                    parent, sign = parents[0]
                    if sign > 0:
                        formula = parent
                    else:
                        formula = f"!{parent}"
                else:
                    # Multiple parents - use OR logic by default
                    parent_terms = []
                    for parent, sign in parents:
                        if sign > 0:
                            parent_terms.append(parent)
                        else:
                            parent_terms.append(f"!{parent}")
                    formula = " | ".join(parent_terms)                
                bnet_lines.append(f"{node}, {formula}")
            else:
                # No parents - input node
                bnet_lines.append(f"{node}, {node}")        
        bnet_content = "\n".join(bnet_lines)        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(bnet_content)        
        return bnet_content
    
    def convert_to_sif(self, learned_edges, output_file=None):
        """Convert learned network to SIF format"""
        
        sif_lines = []
        
        for edge in learned_edges:
            source = edge['source'].replace('TR:', '').replace('DA:', '').replace('DV:', '')
            target = edge['target'].replace('TR:', '').replace('DA:', '').replace('DV:', '')
            
            # Infer sign from correlation or treatment effect
            sign = 1  # Default positive
            if 'correlation' in edge['details']:
                sign = 1 if edge['details']['correlation'] > 0 else -1
            
            sif_lines.append(f"{source}\t{sign}\t{target}")
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write('\n'.join(sif_lines))
        
        return sif_lines
    
    def fit(self, midas_data, prior_pkn_data=None):
        """Main method to learn PKN"""
        
        # Parse data
        data, tr_cols, da_cols, dv_cols = self.parse_midas_data(midas_data)
        
        prior_network = {}
        if prior_pkn_data:
            prior_network = self.parse_pkn_prior(prior_pkn_data)
        
        # Learn network structure
        learned_edges = self.learn_network_structure(data, tr_cols, da_cols, dv_cols, prior_network)
        
        # Store results
        self.learned_edges = learned_edges
        self.data = data
        
        return learned_edges
    
    def get_sif_network(self):
        """Get learned network in SIF format"""
        if not hasattr(self, 'learned_edges'):
            raise ValueError("Must fit model first")
            
        sif_lines = []
        
        for edge in self.learned_edges:
            source = edge['source'].replace('TR:', '').replace('DA:', '').replace('DV:', '').replace('i', '')
            target = edge['target'].replace('TR:', '').replace('DA:', '').replace('DV:', '').replace('i', '')
            
            # Infer sign from statistical evidence
            sign = self.infer_edge_sign(self.data, edge['source'], edge['target'])
            
            sif_lines.append(f"{source}\t{sign}\t{target}")
        
        return sif_lines

# Example usage
if __name__ == "__main__":
    midas_data = "data/ToyModel/ToyModel.csv"
    truth_pkn_file = "data/ToyModel/ToyModel.bnet"
    prior_pkn = "data/ToyModel/90_Modified/ToyModel.bnet"

    # Create and fit model
    model = HybridBayesianPKN(
        alpha=0.1,  # Less stringent for small dataset
        biological_prior_weight=0.9,  # Moderate weight for biological prior
        structure_penalty=0.1
    )
    
    # Learn network
    learned_edges = model.fit(midas_data, prior_pkn)
    
    # Get SIF format
    sif_network = model.get_sif_network()

    model.convert_to_bnet("output/bayesian.bnet")
    
    print("Learned Network (SIF format):")
    print("Source\tSign\tTarget")
    for line in sif_network:
        print(line)
    
    print(f"\nTotal edges learned: {len(sif_network)}")
    
    # Print detailed scores for top edges
    print("\nTop edges with detailed scores:")
    for i, edge in enumerate(model.learned_edges[:10]):
        print(f"{edge['source']} -> {edge['target']}: Score={edge['score']:.3f}")
        print(f"  Statistical: {edge['details']['statistical']:.3f}")
        print(f"  Causal: {edge['details']['causal']:.3f}")
        print(f"  Biological: {edge['details']['biological']:.3f}")
        print(f"  P-value: {edge['details']['p_value']:.3f}")
        print()        
    
    from tools.comparison import AttractorAnalysis
    ori_primes = "data/ToyModel/ToyModel.bnet"
    # compared_bnet = "data/ToyModel/ToyModel.bnet"
    compared_bnet = "output/bayesian.bnet"
    analysis = AttractorAnalysis(ori_primes, compared_bnet)
    results = analysis.comparison()
    print(results)