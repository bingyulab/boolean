import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from scipy import stats
from scipy.stats import pearsonr
import networkx as nx


def save_as_sif(graph: nx.DiGraph, path: str):
    """
    SIF format: each line 'source <interaction> target'
    We'll use '+' for activation, '-' for inhibition.
    """
    with open(path, 'w') as fh:
        for u, v, data in graph.edges(data=True):
            t = data.get('type', data.get('sign', 'activation'))
            if t == 'inhibition' or data.get('sign', 0) == -1:
                fh.write(f"{u}\t-\t{v}\n")
            else:
                fh.write(f"{u}\t+\t{v}\n")
    print(f"Saved SIF to: {path}")

def save_as_bnet(graph: nx.DiGraph, path: str):
    """
    Simple bnet-like CSV: target,factors
    Use '|' to separate OR terms; prefix '!' for inhibition.
    """
    # get all targets that have parents
    targets = sorted({t for _, t in graph.edges()})
    rows = []
    for t in targets:
        parents = []
        for u in sorted(graph.predecessors(t)):
            typ = graph[u][t].get('type', 'activation')
            if typ == 'inhibition' or graph[u][t].get('sign', 0) == -1:
                parents.append('!' + str(u))
            else:
                parents.append(str(u))
        rows.append({'target': t, 'factors': ' | '.join(parents) if parents else ''})
    out = pd.DataFrame(rows)
    out.to_csv(path, index=False)
    print(f"Saved BNET-like CSV to: {path}")


import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from scipy import stats
from scipy.stats import pearsonr
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

class PKNReconstruction:
    def __init__(self, experiment_file, truth_pkn_file, candidate_pkn_file):
        """
        Enhanced PKN Reconstruction Algorithm
        
        Args:
            experiment_file: Path to experimental data CSV
            truth_pkn_file: Path to truth PKN BNET file
            candidate_pkn_file: Path to candidate PKN BNET file
        """
        self.experiment_file = experiment_file
        self.truth_pkn_file = truth_pkn_file
        self.candidate_pkn_file = candidate_pkn_file
        
        # Initialize components
        self.data = None
        self.truth_graph = None
        self.candidate_graph = None
        self.stimuli = []
        self.inhibitors = []
        self.readouts = []
        
        # Results storage
        self.intervention_effects = {}
        self.candidate_edges = {}
        self.final_graph = None
        self.confidence_scores = {}
        
    def load_data(self):
        """Stage 1: Adaptive Data Preprocessing"""
        print("Stage 1: Loading and preprocessing data...")
        
        # Load experimental data
        self.data = pd.read_csv(self.experiment_file)
        
        # Identify column types - handle extra 'i' in inhibitors
        all_tr_cols = [col for col in self.data.columns if col.startswith('TR:')]
        self.inhibitors = []
        self.stimuli = []
        
        for col in all_tr_cols:
            if 'mock:CellLine' in col:
                continue
            # Check if this is an inhibitor (ends with 'i' but not just 'i')
            if col.endswith('i') and len(col) > 4:  # More than just 'TR:i'
                # Check if the base name exists (e.g., TR:mek12 for TR:mek12i)
                base_name = col[:-1]  # Remove the 'i'
                if base_name not in all_tr_cols:
                    # This is an inhibitor
                    self.inhibitors.append(col)
                else:
                    # This might be a stimulus that naturally ends in 'i'
                    self.stimuli.append(col)
            else:
                self.stimuli.append(col)
        
        self.readouts = [col for col in self.data.columns if col.startswith('DV:')]
        
        print(f"Found {len(self.stimuli)} stimuli: {self.stimuli}")
        print(f"Found {len(self.inhibitors)} inhibitors: {self.inhibitors}")
        print(f"Found {len(self.readouts)} readouts: {self.readouts}")
        
        # Keep only intervention rows (DA time != 0)
        da_cols = [col for col in self.data.columns if col.startswith('DA:')]
        if da_cols:
            intervention_mask = self.data[da_cols].sum(axis=1) > 0
            self.data = self.data[intervention_mask].reset_index(drop=True)
            print(f"Kept {len(self.data)} intervention rows")
        
        # Handle NaN values
        for col in self.readouts:
            self.data[col] = self.data[col].fillna(0)
        
        # Quality assessment
        self._assess_data_quality()
        
    def _assess_data_quality(self):
        """Assess data quality per node"""
        self.data_quality = {}
        for readout in self.readouts:
            values = self.data[readout].dropna()
            if len(values) > 0:
                self.data_quality[readout] = {
                    'completeness': len(values) / len(self.data),
                    'variance': values.var(),
                    'range': values.max() - values.min(),
                    'outlier_fraction': len(values[(np.abs(stats.zscore(values)) > 3)]) / len(values)
                }
            else:
                self.data_quality[readout] = {'completeness': 0, 'variance': 0, 'range': 0, 'outlier_fraction': 0}
    
    def load_pkn_networks(self):
        """Load PKN networks from BNET files"""
        print("Loading PKN networks...")
        
        self.truth_graph = self._parse_bnet(self.truth_pkn_file)
        self.candidate_graph = self._parse_bnet(self.candidate_pkn_file)
        
        print(f"Truth PKN: {self.truth_graph.number_of_nodes()} nodes, {self.truth_graph.number_of_edges()} edges")
        print(f"Candidate PKN: {self.candidate_graph.number_of_nodes()} nodes, {self.candidate_graph.number_of_edges()} edges")
    
    def _parse_bnet(self, filename):
        """Parse BNET file format"""
        graph = nx.DiGraph()
        
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()[1:]
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if ',' in line:
                    parts = line.split(',')
                    target = parts[0].strip()
                    factors = parts[1].strip()
                    
                    # Parse factors (assuming OR logic with | and inhibition with !)
                    if factors:
                        # Split by | for OR relationships
                        or_factors = factors.split('|')
                        for factor in or_factors:
                            factor = factor.strip()
                            if factor.startswith('!'):
                                source = factor[1:].strip()
                                edge_type = 'inhibition'
                            else:
                                source = factor.strip()
                                edge_type = 'activation'
                            
                            if source and target:
                                graph.add_edge(source, target, type=edge_type)
        except FileNotFoundError:
            print(f"Warning: {filename} not found. Creating empty graph.")
        
        return graph
    
    def compute_intervention_effects(self):
        """Stage 2: Ensemble Effect Discovery"""
        print("Stage 2: Computing intervention effects...")
        
        self.intervention_effects = {}
        
        # For each stimulus-readout pair
        for stimulus in self.stimuli:
            for readout in self.readouts:
                # Clean column names
                stim_name = stimulus.replace('TR:', '')
                read_name = readout.replace('DV:', '')
                
                effects = self._compute_multi_method_effects(stimulus, readout)
                self.intervention_effects[(stim_name, read_name)] = effects
        
        # Compute mediator effects
        self._compute_mediator_effects()
    
    def _compute_multi_method_effects(self, stimulus, readout):
        """Compute intervention effects using multiple methods"""
        effects = {}
        
        # Method 1: Direct effect (controlling for other stimuli)
        stim_on = self.data[self.data[stimulus] == 1]
        stim_off = self.data[self.data[stimulus] == 0]
        
        if len(stim_on) > 0 and len(stim_off) > 0:
            mean_on = stim_on[readout].mean()
            mean_off = stim_off[readout].mean()
            effects['direct'] = mean_on - mean_off
            
            # Statistical significance
            if len(stim_on) > 1 and len(stim_off) > 1:
                _, p_val = stats.ttest_ind(stim_on[readout], stim_off[readout])
                effects['p_value'] = p_val
            else:
                effects['p_value'] = 1.0
        else:
            effects['direct'] = 0
            effects['p_value'] = 1.0
        
        # Method 2: Correlation-based
        if len(self.data) > 2:
            corr, _ = pearsonr(self.data[stimulus], self.data[readout])
            effects['correlation'] = corr if not np.isnan(corr) else 0
        else:
            effects['correlation'] = 0
        
        # Method 3: Mutual information (simplified)
        effects['mutual_info'] = self._compute_mutual_info(self.data[stimulus], self.data[readout])
        
        return effects
    
    def _compute_mutual_info(self, x, y):
        """Simplified mutual information calculation"""
        try:
            # Discretize continuous variables
            x_disc = np.digitize(x, np.percentile(x, [25, 50, 75]))
            y_disc = np.digitize(y, np.percentile(y, [25, 50, 75]))
            
            # Compute mutual information
            xy = np.vstack([x_disc, y_disc]).T
            unique_xy, counts_xy = np.unique(xy, axis=0, return_counts=True)
            p_xy = counts_xy / len(xy)
            
            unique_x, counts_x = np.unique(x_disc, return_counts=True)
            p_x = counts_x / len(x_disc)
            
            unique_y, counts_y = np.unique(y_disc, return_counts=True)
            p_y = counts_y / len(y_disc)
            
            mi = 0
            for i, (x_val, y_val) in enumerate(unique_xy):
                p_joint = p_xy[i]
                p_x_marg = p_x[x_val]
                p_y_marg = p_y[y_val]
                if p_joint > 0 and p_x_marg > 0 and p_y_marg > 0:
                    mi += p_joint * np.log2(p_joint / (p_x_marg * p_y_marg))
            
            return mi
        except:
            return 0
    
    def _compute_mediator_effects(self):
        """Test for mediator relationships"""
        print("Computing mediator effects...")
        
        # For each potential mediator relationship
        for (stim, readout), effects in self.intervention_effects.items():
            if effects['direct'] != 0:  # Only for non-zero effects
                # Test potential mediators (other readouts)
                mediator_scores = {}
                
                stim_col = f'TR:{stim}'
                readout_col = f'DV:{readout}'
                
                if stim_col in self.data.columns and readout_col in self.data.columns:
                    for potential_mediator in self.readouts:
                        if potential_mediator != readout_col:
                            mediator_score = self._test_mediation(stim_col, potential_mediator, readout_col)
                            med_name = potential_mediator.replace('DV:', '')
                            mediator_scores[med_name] = mediator_score
                
                effects['mediators'] = mediator_scores
    
    def _test_mediation(self, stimulus, mediator, outcome):
        """Simple mediation test"""
        try:
            # Direct path coefficient (X -> Y)
            x = self.data[stimulus]
            y = self.data[outcome]
            m = self.data[mediator]
            
            # Simple regression coefficients
            direct_coef = np.corrcoef(x, y)[0, 1] if len(set(x)) > 1 else 0
            
            # Mediated path (X -> M, M -> Y controlling for X)
            if len(set(x)) > 1 and len(set(m)) > 1:
                x_to_m = np.corrcoef(x, m)[0, 1]
                
                # Partial correlation M -> Y controlling for X
                if len(set(y)) > 1:
                    # Simplified partial correlation
                    xy_corr = np.corrcoef(x, y)[0, 1]
                    xm_corr = np.corrcoef(x, m)[0, 1]
                    my_corr = np.corrcoef(m, y)[0, 1]
                    
                    if not (np.isnan(xy_corr) or np.isnan(xm_corr) or np.isnan(my_corr)):
                        partial_corr = (my_corr - xy_corr * xm_corr) / (np.sqrt(1 - xy_corr**2) * np.sqrt(1 - xm_corr**2))
                        mediation_strength = abs(x_to_m * partial_corr)
                        return mediation_strength
            
            return 0
        except:
            return 0
    
    def adaptive_parent_selection(self):
        """Stage 3: Adaptive Parent Selection using regularized regression"""
        print("Stage 3: Adaptive parent selection...")
        
        self.parent_candidates = {}
        
        for readout in self.readouts:
            read_name = readout.replace('DV:', '')
            print(f"  Selecting parents for {read_name}")
            
            # Prepare features (stimuli and other readouts)
            feature_cols = self.stimuli + [r for r in self.readouts if r != readout]
            
            # Get data with no missing values for this readout
            valid_mask = self.data[readout].notna()
            if valid_mask.sum() < 5:  # Need minimum samples
                self.parent_candidates[read_name] = []
                continue
            
            X = self.data[valid_mask][feature_cols].fillna(0)
            y = self.data[valid_mask][readout]
            
            if len(set(y)) > 1:  # Need variability in target
                parents = self._robust_parent_selection(X, y, read_name)
                self.parent_candidates[read_name] = parents
            else:
                self.parent_candidates[read_name] = []
    
    def _robust_parent_selection(self, X, y, target_name):
        """Robust parent selection with multiple methods"""
        n_features = X.shape[1]
        if n_features == 0:
            return []
        
        parent_scores = {}
        
        # Method 1: LASSO with CV
        try:
            lasso = LassoCV(cv=min(5, len(X)//2), random_state=42, max_iter=1000)
            lasso.fit(X, y)
            
            for i, feature in enumerate(X.columns):
                if abs(lasso.coef_[i]) > 1e-6:
                    parent_scores[feature] = abs(lasso.coef_[i])
        except:
            pass
        
        # Method 2: Correlation-based selection
        for feature in X.columns:
            if len(set(X[feature])) > 1:
                corr = abs(np.corrcoef(X[feature], y)[0, 1])
                if not np.isnan(corr):
                    if feature in parent_scores:
                        parent_scores[feature] = max(parent_scores[feature], corr)
                    else:
                        parent_scores[feature] = corr
        
        # Method 3: PKN-informed scoring (boost edges in candidate PKN)
        if self.candidate_graph:
            target_clean = target_name
            for feature in X.columns:
                feature_clean = feature.replace('TR:', '').replace('DV:', '')
                if self.candidate_graph.has_edge(feature_clean, target_clean):
                    if feature in parent_scores:
                        parent_scores[feature] *= 1.5  # Boost PKN edges
                    else:
                        parent_scores[feature] = 0.1  # Give minimal score to PKN edges
        
        # Select top parents
        if parent_scores:
            # Dynamic threshold based on data quality
            quality = self.data_quality.get(f'DV:{target_name}', {})
            base_threshold = 0.1
            quality_factor = quality.get('completeness', 0.5) * quality.get('variance', 0.5)
            threshold = base_threshold * quality_factor
            
            selected_parents = [p for p, score in parent_scores.items() if score > threshold]
            return selected_parents[:5]  # Limit to top 5 parents
        
        return []

    def hybrid_causal_structure_learning(self):
        """Stage 4: Hybrid Causal Structure Learning with Inhibitor Detection"""
        print("Stage 4: Hybrid causal structure learning...")
        
        self.candidate_edges = {}
        
        # Method 1: Intervention-based edge scoring (including inhibitors)
        for (stim, readout), effects in self.intervention_effects.items():
            edge_key = (stim, readout)
            
            # Combine multiple effect measures
            direct_effect = effects.get('direct', 0)
            correlation = effects.get('correlation', 0)
            mutual_info = effects.get('mutual_info', 0)
            p_value = effects.get('p_value', 1.0)
            
            # Ensemble scoring
            ensemble_score = (
                0.4 * abs(direct_effect) + 
                0.3 * abs(correlation) + 
                0.3 * mutual_info
            ) * (1 - p_value)  # Weight by significance
            
            if ensemble_score > 0.1:  # Threshold for candidate edges
                self.candidate_edges[edge_key] = {
                    'score': ensemble_score,
                    'type': 'activation' if direct_effect > 0 else 'inhibition',
                    'evidence': effects
                }
        
        # Method 2: Inhibitor analysis
        self._analyze_inhibitor_effects()
        
        # Method 3: Parent selection based edges
        for readout, parents in self.parent_candidates.items():
            for parent in parents:
                parent_clean = parent.replace('TR:', '').replace('DV:', '')
                edge_key = (parent_clean, readout)
                
                if edge_key not in self.candidate_edges:
                    self.candidate_edges[edge_key] = {
                        'score': 0.2,  # Base score for parent-selected edges
                        'type': 'activation',
                        'evidence': {'method': 'parent_selection'}
                    }
                else:
                    # Boost score if found by multiple methods
                    self.candidate_edges[edge_key]['score'] *= 1.3
    
    def _analyze_inhibitor_effects(self):
        """Analyze inhibitor effects to identify regulated pathways"""
        print("  Analyzing inhibitor effects...")
        
        for inhibitor in self.inhibitors:
            # Extract target name (remove 'i' and 'TR:')
            target_name = inhibitor.replace('TR:', '').replace('i', '')
            
            # Find experiments with/without this inhibitor
            inhib_on = self.data[self.data[inhibitor] == 1]
            inhib_off = self.data[self.data[inhibitor] == 0]
            
            if len(inhib_on) > 0 and len(inhib_off) > 0:
                # Compare readouts when inhibitor is on vs off
                for readout in self.readouts:
                    readout_name = readout.replace('DV:', '')
                    
                    mean_with_inhib = inhib_on[readout].mean()
                    mean_without_inhib = inhib_off[readout].mean()
                    
                    inhibitor_effect = mean_without_inhib - mean_with_inhib  # Positive means inhibitor reduces signal
                    
                    if abs(inhibitor_effect) > 0.1:  # Significant effect
                        # This suggests target_name regulates readout_name
                        edge_key = (target_name, readout_name)
                        
                        # Compute significance
                        if len(inhib_on) > 1 and len(inhib_off) > 1:
                            _, p_val = stats.ttest_ind(inhib_on[readout].dropna(), 
                                                     inhib_off[readout].dropna())
                        else:
                            p_val = 0.5
                        
                        score = abs(inhibitor_effect) * (1 - p_val)
                        
                        if edge_key in self.candidate_edges:
                            # Boost existing edge score
                            self.candidate_edges[edge_key]['score'] = max(
                                self.candidate_edges[edge_key]['score'], score
                            )
                            self.candidate_edges[edge_key]['evidence']['inhibitor_effect'] = inhibitor_effect
                        else:
                            # Add new edge from inhibitor analysis
                            self.candidate_edges[edge_key] = {
                                'score': score,
                                'type': 'activation' if inhibitor_effect > 0 else 'inhibition',
                                'evidence': {
                                    'method': 'inhibitor_analysis',
                                    'inhibitor_effect': inhibitor_effect,
                                    'inhibitor': inhibitor,
                                    'p_value': p_val
                                }
                            }
                        
                        print(f"    {target_name} -> {readout_name}: effect={inhibitor_effect:.3f}, p={p_val:.3f}")
    
    def detect_intermediate_nodes(self):
        """Stage 4.5: Detect intermediate nodes in pathways"""
        print("Stage 4.5: Detecting intermediate nodes...")
        
        self.intermediate_pathways = {}
        
        # For each direct connection found, test for intermediate nodes
        direct_edges = list(self.candidate_edges.keys())
        
        for (source, target) in direct_edges:
            source_col = f'TR:{source}' if f'TR:{source}' in self.data.columns else f'DV:{source}'
            target_col = f'DV:{target}'
            
            if source_col in self.data.columns and target_col in self.data.columns:
                # Test all potential intermediates
                potential_intermediates = []
                
                for readout in self.readouts:
                    intermediate = readout.replace('DV:', '')
                    if intermediate != source and intermediate != target:
                        mediation_score = self._test_pathway_mediation(source_col, readout, target_col)
                        
                        if mediation_score > 0.3:  # Significant mediation
                            potential_intermediates.append({
                                'node': intermediate,
                                'mediation_score': mediation_score,
                                'pathway': f"{source} -> {intermediate} -> {target}"
                            })
                
                if potential_intermediates:
                    # Sort by mediation score
                    potential_intermediates.sort(key=lambda x: x['mediation_score'], reverse=True)
                    self.intermediate_pathways[(source, target)] = potential_intermediates[:3]  # Top 3
                    
                    print(f"  Found intermediates for {source} -> {target}:")
                    for intermediate in potential_intermediates[:2]:
                        print(f"    {intermediate['pathway']} (score: {intermediate['mediation_score']:.3f})")
        
    # ------------ NEW: roles & candidate graph ------------
    def _derive_roles(self):
        """Define sources (stimuli), targets (readouts), and allowed middle nodes."""
        self.sources = [c.replace('TR:', '') if c.startswith('TR:') else c for c in self.stimuli]
        self.targets = [c.replace('DV:', '') if c.startswith('DV:') else c for c in self.readouts]
        # middle nodes are readouts that can be used as intermediates
        self.middle_nodes = list(sorted(set(self.targets)))  # keep it simple: DVs as intermediates

    def _build_candidate_graph(self, min_edge_score=0.10):
        """
        Build a directed graph of candidate edges using self.candidate_edges.
        Keeps edges with score >= min_edge_score.
        """
        G = nx.DiGraph()
        for (u, v), info in self.candidate_edges.items():
            if info.get('score', 0.0) >= min_edge_score:
                G.add_node(u); G.add_node(v)
                G.add_edge(u, v, score=float(info['score']),
                        type=info.get('type', 'activation'))
        return G

    # ------------ NEW: path enumeration ------------
    def _enumerate_sr_paths(self, G, s, r, max_len=3, top_k=25):
        """
        Enumerate simple paths s->r with length <= max_len through allowed middle nodes.
        Score of a path = geometric mean of edge scores.
        Returns list of dicts: {'nodes': [...], 'edges': [(u,v),...], 'score': float}
        """
        if s not in G or r not in G:
            return []
        results = []

        def dfs(current, target, path_nodes, path_edges, depth):
            if depth > max_len:
                return
            if current == target:
                if path_edges:
                    # geometric mean of edge scores
                    es = [G[u][v]['score'] for (u, v) in path_edges]
                    gm = float(np.exp(np.mean(np.log(np.maximum(es, 1e-8)))))
                    results.append({'nodes': path_nodes[:], 'edges': path_edges[:], 'score': gm})
                return
            # expand
            for nbr in G.successors(current):
                # allow only middle nodes (DVs) as intermediates, except the final target
                if nbr != target and nbr not in self.middle_nodes:
                    continue
                if nbr in path_nodes:
                    continue  # avoid cycles
                path_nodes.append(nbr)
                path_edges.append((current, nbr))
                dfs(nbr, target, path_nodes, path_edges, depth + 1)
                path_edges.pop()
                path_nodes.pop()

        dfs(s, r, [s], [], 0)
        # keep top-k by path score
        results.sort(key=lambda d: d['score'], reverse=True)
        return results[:top_k]

    # ------------ NEW: ILP selection of path cover ------------
    def _select_paths_via_ilp(self, all_paths, lambda_edges=0.15, max_paths_per_pair=2):
        """
        Solve an ILP to select a sparse subset of paths.
        Variables:
        y_p ∈ {0,1} for each path p
        x_e ∈ {0,1} for each edge e  (edge chosen by any selected path)
        Objective:
        maximize  sum_p w_p * y_p  -  lambda_edges * sum_e x_e
        with constraints: y_p <= x_e for all edges e in path p
                            at most max_paths_per_pair per (s,r)
        Returns: selected_edges (set of (u,v)), selected_paths (list of path dicts)
        """
        from pulp import LpProblem, LpVariable, LpMaximize, LpBinary, lpSum, PULP_CBC_CMD

        # collect unique edges, index paths
        edges = set()
        for p in all_paths:
            for e in p['edges']:
                edges.add(e)
        edges = sorted(edges)
        edge_idx = {e: i for i, e in enumerate(edges)}

        # group paths by (s,r)
        pair_to_paths = {}
        for i, p in enumerate(all_paths):
            s, r = p['nodes'][0], p['nodes'][-1]
            pair_to_paths.setdefault((s, r), []).append((i, p))

        # define ILP
        prob = LpProblem("StimulusReadoutPathCover", LpMaximize)
        y = [LpVariable(f"y_{i}", lowBound=0, upBound=1, cat=LpBinary) for i in range(len(all_paths))]
        x = [LpVariable(f"x_{i}", lowBound=0, upBound=1, cat=LpBinary) for i in range(len(edges))]

        # objective weights: path score * demand weight
        # demand weight from direct effect (normalized to [0,1])
        def _effect_weight(s, r):
            eff = abs(self.intervention_effects.get((s, r), {}).get('direct', 0.0))
            return float(min(1.0, eff / (eff + 1e-8)))

        obj_terms = []
        for i, p in enumerate(all_paths):
            s, r = p['nodes'][0], p['nodes'][-1]
            w = p['score'] * _effect_weight(s, r)
            obj_terms.append(w * y[i])
        # penalty on edges
        edge_penalty = lambda_edges * lpSum(x)

        prob += lpSum(obj_terms) - edge_penalty

        # linking constraints: y_p <= x_e for all edges in path p
        for i, p in enumerate(all_paths):
            for e in p['edges']:
                prob += y[i] <= x[edge_idx[e]]

        # per-pair cap: at most max_paths_per_pair
        for pair, plist in pair_to_paths.items():
            prob += lpSum([y[i] for i, _ in plist]) <= max_paths_per_pair

        # solve
        prob.solve(PULP_CBC_CMD(msg=False))

        # collect solution
        selected_paths = []
        for i, p in enumerate(all_paths):
            if y[i].value() is not None and y[i].value() > 0.5:
                selected_paths.append(p)

        selected_edges = set()
        for e, var in zip(edges, x):
            if var.value() is not None and var.value() > 0.5:
                selected_edges.add(e)

        return selected_edges, selected_paths

    # ------------ NEW: public runner ------------
    def build_pkn_via_path_cover(self,
                                min_edge_score=0.10,
                                max_len=3,
                                top_k_per_pair=25,
                                lambda_edges=0.15,
                                max_paths_per_pair=2,
                                overwrite_final=True):
        """
        Construct a PKN by enumerating s->r paths and selecting a sparse subset via ILP.
        Produces self.final_graph and updates confidence/type from candidate edges.
        """
        print("Path-cover PKN construction:")
        # Make sure candidate_edges exist
        if not self.candidate_edges:
            raise RuntimeError("candidate_edges is empty. Run compute_intervention_effects(), "
                            "adaptive_parent_selection(), and hybrid_causal_structure_learning() first.")

        self._derive_roles()
        Gc = self._build_candidate_graph(min_edge_score=min_edge_score)

        # enumerate paths for all (s, r) that have non-trivial effect
        all_paths = []
        sr_pairs = []
        for s in self.sources:
            for r in self.targets:
                eff = abs(self.intervention_effects.get((s, r), {}).get('direct', 0.0))
                if eff < 1e-6:
                    continue
                paths = self._enumerate_sr_paths(Gc, s, r, max_len=max_len, top_k=top_k_per_pair)
                if paths:
                    all_paths.extend(paths)
                    sr_pairs.append((s, r, len(paths)))
        print(f"  Enumerated {len(all_paths)} paths across {len(sr_pairs)} (stimulus, readout) pairs.")

        # select paths via ILP
        selected_edges, selected_paths = self._select_paths_via_ilp(
            all_paths,
            lambda_edges=lambda_edges,
            max_paths_per_pair=max_paths_per_pair
        )
        print(f"  Selected {len(selected_paths)} paths using {len(selected_edges)} edges.")

        # build final_graph
        if overwrite_final or self.final_graph is None:
            self.final_graph = nx.DiGraph()

        for (u, v) in selected_edges:
            # decide sign from candidate_edges (majority across methods)
            typ = 'activation'
            if (u, v) in self.candidate_edges:
                typ = self.candidate_edges[(u, v)].get('type', 'activation')
                conf = self.candidate_edges[(u, v)].get('score', 0.2)
            else:
                # fallback: positive by default
                conf = 0.2
            self.final_graph.add_node(u); self.final_graph.add_node(v)
            self.final_graph.add_edge(u, v, type=typ, confidence=conf)

        print(f"Final (path-cover) graph: {self.final_graph.number_of_nodes()} nodes, {self.final_graph.number_of_edges()} edges.")

    def _test_pathway_mediation(self, source_col, mediator_col, target_col):
        """Enhanced mediation test for pathway detection"""
        try:
            X = self.data[source_col].values
            M = self.data[mediator_col].values  
            Y = self.data[target_col].values
            
            # Remove rows with NaN
            valid_mask = ~(np.isnan(X) | np.isnan(M) | np.isnan(Y))
            X, M, Y = X[valid_mask], M[valid_mask], Y[valid_mask]
            
            if len(X) < 5:
                return 0
            
            # Sobel test approximation
            # Path a: X -> M
            if len(set(X)) > 1 and len(set(M)) > 1:
                a = np.corrcoef(X, M)[0, 1]
            else:
                a = 0
                
            # Path b: M -> Y (controlling for X)
            if len(set(M)) > 1 and len(set(Y)) > 1:
                # Simple partial correlation
                rm_y = np.corrcoef(M, Y)[0, 1] if not np.isnan(np.corrcoef(M, Y)[0, 1]) else 0
                rx_y = np.corrcoef(X, Y)[0, 1] if not np.isnan(np.corrcoef(X, Y)[0, 1]) else 0
                rx_m = np.corrcoef(X, M)[0, 1] if not np.isnan(np.corrcoef(X, M)[0, 1]) else 0
                
                denominator = np.sqrt((1 - rx_y**2) * (1 - rx_m**2))
                if denominator > 1e-6:
                    b = (rm_y - rx_y * rx_m) / denominator
                else:
                    b = 0
            else:
                b = 0
            
            # Mediation effect
            mediation_effect = abs(a * b)
            
            # Direct effect
            direct_effect = abs(rx_y - a * b) if 'rx_y' in locals() else 0
            
            # Mediation strength (proportion mediated)
            if direct_effect + mediation_effect > 1e-6:
                mediation_strength = mediation_effect / (direct_effect + mediation_effect)
                return mediation_strength
            else:
                return mediation_effect
                
        except Exception as e:
            return 0
    
    def dynamic_pkn_refinement(self):
        """Stage 5: Dynamic PKN Refinement"""
        print("Stage 5: Dynamic PKN refinement...")
        
        # Initialize with candidate PKN
        self.final_graph = self.candidate_graph.copy() if self.candidate_graph else nx.DiGraph()
        
        # Add high-confidence edges from data
        for (source, target), edge_info in self.candidate_edges.items():
            if edge_info['score'] > 0.3:  # High confidence threshold
                edge_type = edge_info['type']
                self.final_graph.add_edge(source, target, 
                                        type=edge_type, 
                                        confidence=edge_info['score'],
                                        source='data')
        
        # Remove low-evidence edges from PKN
        edges_to_remove = []
        for source, target in self.final_graph.edges():
            edge_key = (source, target)
            if edge_key not in self.candidate_edges:
                # Check if this PKN edge has any data support
                data_support = self._check_data_support(source, target)
                if data_support < 0.05:  # Very low support
                    edges_to_remove.append((source, target))
        
        for edge in edges_to_remove:
            self.final_graph.remove_edge(*edge)
            print(f"Removed low-support PKN edge: {edge}")
        
        print(f"Final graph: {self.final_graph.number_of_nodes()} nodes, {self.final_graph.number_of_edges()} edges")
    
    def _check_data_support(self, source, target):
        """Check if a PKN edge has support in the data"""
        # Look for evidence in intervention effects
        edge_key = (source, target)
        if edge_key in self.intervention_effects:
            effects = self.intervention_effects[edge_key]
            return abs(effects.get('direct', 0)) + abs(effects.get('correlation', 0))
        return 0
    
    def validate_and_score(self):
        """Stage 6: Validation and confidence scoring"""
        print("Stage 6: Validation and confidence scoring...")
        
        self.validation_results = {}
        
        # Compute confidence scores for each edge
        for source, target in self.final_graph.edges():
            confidence = self._compute_edge_confidence(source, target)
            self.confidence_scores[(source, target)] = confidence
            self.final_graph[source][target]['confidence'] = confidence
        
        # Overall network metrics
        self.validation_results['total_edges'] = self.final_graph.number_of_edges()
        self.validation_results['avg_confidence'] = np.mean(list(self.confidence_scores.values()))
        
        # Compare with truth if available
        if self.truth_graph:
            self._compare_with_truth()
    
    def _compute_edge_confidence(self, source, target):
        """Compute confidence score for an edge"""
        confidence = 0.5  # Base confidence
        
        # Check data evidence
        edge_key = (source, target)
        if edge_key in self.candidate_edges:
            confidence = self.candidate_edges[edge_key]['score']
        elif edge_key in self.intervention_effects:
            effects = self.intervention_effects[edge_key]
            confidence = abs(effects.get('direct', 0)) * (1 - effects.get('p_value', 1))
        
        # PKN consistency bonus
        # if self.candidate_graph and self.candidate_graph.has_edge(source, target):
        #     confidence *= 1.2
        
        return min(confidence, 1.0)
    
    def _compare_with_truth(self):
        """Compare reconstructed network with ground truth"""
        if not self.truth_graph:
            return
        
        true_edges = set(self.truth_graph.edges())
        pred_edges = set(self.final_graph.edges())
        
        tp = len(true_edges.intersection(pred_edges))
        fp = len(pred_edges - true_edges)
        fn = len(true_edges - pred_edges)
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        self.validation_results.update({
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'true_positives': tp,
            'false_positives': fp,
            'false_negatives': fn
        })
        
        print(f"Validation results:")
        print(f"  Precision: {precision:.3f}")
        print(f"  Recall: {recall:.3f}")
        print(f"  F1 Score: {f1:.3f}")
    
    def generate_report(self):
        """Generate comprehensive analysis report"""
        print("\n" + "="*60)
        print("PKN RECONSTRUCTION ANALYSIS REPORT")
        print("="*60)
        
        print(f"\nDATA SUMMARY:")
        print(f"  Total experiments: {len(self.data)}")
        print(f"  Stimuli: {len(self.stimuli)}")
        print(f"  Inhibitors: {len(self.inhibitors)}")
        print(f"  Readouts: {len(self.readouts)}")
        
        print(f"\nNETWORK COMPARISON:")
        if self.candidate_graph:
            print(f"  Candidate PKN edges: {self.candidate_graph.number_of_edges()}")
        if self.truth_graph:
            print(f"  Truth PKN edges: {self.truth_graph.number_of_edges()}")
        print(f"  Final network edges: {self.final_graph.number_of_edges()}")
        
        print(f"\nVALIDATION METRICS:")
        if 'precision' in self.validation_results:
            print(f"  Precision: {self.validation_results['precision']:.3f}")
            print(f"  Recall: {self.validation_results['recall']:.3f}")
            print(f"  F1 Score: {self.validation_results['f1_score']:.3f}")
        
        print(f"\nTOP CONFIDENT EDGES:")
        sorted_edges = sorted(self.confidence_scores.items(), 
                            key=lambda x: x[1], reverse=True)[:10]
        for (source, target), conf in sorted_edges:
            print(f"  {source} -> {target}: {conf:.3f}")
        
        print(f"\nINTERMEDIATE PATHWAYS DETECTED:")
        if hasattr(self, 'intermediate_pathways'):
            for (source, target), intermediates in self.intermediate_pathways.items():
                print(f"  Direct: {source} -> {target}")
                for intermediate in intermediates[:2]:  # Show top 2
                    print(f"    Via: {intermediate['pathway']} (score: {intermediate['mediation_score']:.3f})")
        
        print(f"\nINHIBITOR EFFECTS:")
        inhibitor_edges = [(k, v) for k, v in self.candidate_edges.items() 
                          if v.get('evidence', {}).get('method') == 'inhibitor_analysis']
        for (source, target), info in inhibitor_edges[:5]:
            effect = info['evidence'].get('inhibitor_effect', 0)
            print(f"  {source} -> {target}: effect={effect:.3f}")
        
        print(f"\nCAUSAL DISCOVERY RESULTS:")
        if hasattr(self, 'pc_edges'):
            pc_only_edges = [(k, v) for k, v in self.pc_edges.items() 
                            if k not in self.intervention_effects]
            print(f"  PC algorithm found {len(pc_only_edges)} additional edges")
            for (source, target), conf in list(pc_only_edges)[:3]:
                print(f"    {source} -> {target}: {conf:.3f}")
        
        print(f"\nNOVEL DISCOVERIES:")
        novel_edges = []
        for source, target in self.final_graph.edges():
            if not (self.candidate_graph and self.candidate_graph.has_edge(source, target)):
                conf = self.confidence_scores.get((source, target), 0)
                novel_edges.append(((source, target), conf))
        
        novel_edges.sort(key=lambda x: x[1], reverse=True)
        for (source, target), conf in novel_edges[:5]:
            print(f"  {source} -> {target}: {conf:.3f}")
        
        return {
            'final_graph': self.final_graph,
            'confidence_scores': self.confidence_scores,
            'validation_results': self.validation_results,
            'intervention_effects': self.intervention_effects,
            'intermediate_pathways': getattr(self, 'intermediate_pathways', {}),
            'pc_edges': getattr(self, 'pc_edges', {})
        }    
    
    def run_full_analysis(self):
        """Run the complete PKN reconstruction pipeline"""
        print("Starting Enhanced PKN Reconstruction Analysis...")
        
        # Execute all stages
        self.load_data()
        self.load_pkn_networks()
        self.compute_intervention_effects()
        self.adaptive_parent_selection()
        self.hybrid_causal_structure_learning()
        self.build_pkn_via_path_cover(
                min_edge_score=0.10,   # keep modest candidates
                max_len=3,             # typical is 2–3 to control combinatorics
                top_k_per_pair=20,     # prune per-pair
                lambda_edges=0.15,     # sparsity penalty
                max_paths_per_pair=2   # avoid overfitting per (s,r)
            )
        self.detect_intermediate_nodes()  # New stage
        self.dynamic_pkn_refinement()
        self.validate_and_score()
        
        # Generate final report
        results = self.generate_report()

        save_as_bnet(results['final_graph'], "output/flow.bnet")
        return results


if __name__ == "__main__":
    reconstructor = PKNReconstruction(
        experiment_file="data/ToyModel/ToyModel.csv",
        truth_pkn_file="data/ToyModel/ToyModel.bnet", 
        candidate_pkn_file="data/ToyModel/90_Modified/ToyModel.bnet"
    )

    # Run full analysis
    results = reconstructor.run_full_analysis()

    print(f"results: {results}")
    
    from tools.comparison import AttractorAnalysis
    ori_primes = "data/ToyModel/ToyModel.bnet"
    # compared_bnet = "data/ToyModel/ToyModel.bnet"
    compared_bnet = "output/OwnMethod.bnet"
    analysis = AttractorAnalysis(ori_primes, compared_bnet)
    results = analysis.comparison()
    print(results)