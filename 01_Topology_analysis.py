import networkx as nx
import pandas as pd
import numpy as np
from cellnopt.core import CNOGraph
from cellnopt.boolean import BooleanModel
import caspo
from caspo.core import Graph, Setup
from caspo.learn import learn_ensemble

def load_sif(filename):
    G = nx.DiGraph()
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                source, interaction, target = parts
                G.add_edge(source, target, interaction=interaction)
    return G

# Method 1: Using CellNOpt for topology analysis
def analyze_with_cellnopt(pkn_graph):
    """Analyze network topology using CellNOpt"""
    print("=== CellNOpt Topology Analysis ===")
    
    # Convert NetworkX graph to CellNOpt format
    # Create SIF (Simple Interaction Format) representation
    sif_data = []
    for u, v, data in pkn_graph.edges(data=True):
        interaction = data.get('interaction', 'unknown')
        sif_data.append([u, interaction, v])
    
    sif_df = pd.DataFrame(sif_data, columns=['source', 'interaction', 'target'])
    
    # Basic network statistics
    print(f"Network nodes: {pkn_graph.number_of_nodes()}")
    print(f"Network edges: {pkn_graph.number_of_edges()}")
    print(f"Network density: {nx.density(pkn_graph):.3f}")
    
    # Find key network properties
    try:
        # Strongly connected components
        sccs = list(nx.strongly_connected_components(pkn_graph))
        print(f"Strongly connected components: {len(sccs)}")
        
        # Weakly connected components
        wccs = list(nx.weakly_connected_components(pkn_graph))
        print(f"Weakly connected components: {len(wccs)}")
        
        # Find paths between input/output nodes
        inputs = [node for node in pkn_graph.nodes() if pkn_graph.in_degree(node) == 0]
        outputs = [node for node in pkn_graph.nodes() if pkn_graph.out_degree(node) == 0]
        
        print(f"Input nodes: {inputs}")
        print(f"Output nodes: {outputs}")
        
    except Exception as e:
        print(f"Error in topology analysis: {e}")
    
    return sif_df

# Method 2: Using CASPO for logic-based analysis
def analyze_with_caspo(pkn_graph):
    """Analyze network using CASPO logic-based approach"""
    print("\n=== CASPO Logic-Based Analysis ===")
    
    try:
        # Convert to CASPO format
        nodes = list(pkn_graph.nodes())
        edges = [(u, v) for u, v in pkn_graph.edges()]
        
        # Create CASPO graph
        caspo_graph = Graph.read_sif_lines([f"{u}\t1\t{v}" for u, v in edges])
        
        # Generate Boolean rules based on network topology
        print("Generating Boolean logic rules...")
        
        # Simple rule generation based on network structure
        rules = {}
        for node in pkn_graph.nodes():
            predecessors = list(pkn_graph.predecessors(node))
            if predecessors:
                if len(predecessors) == 1:
                    rules[node] = f"{node} = {predecessors[0]}"
                else:
                    # Use OR logic for multiple inputs (can be customized)
                    rule = f"{node} = " + " | ".join(predecessors)
                    rules[node] = rule
            else:
                rules[node] = f"{node} = {node}"  # Input node
        
        print("Generated Boolean rules:")
        for node, rule in rules.items():
            print(f"  {rule}")
            
        return rules
        
    except Exception as e:
        print(f"Error in CASPO analysis: {e}")
        return {}

# Method 3: Network-based optimization without data
def optimize_network_structure(pkn_graph):
    """Optimize network structure using graph theory metrics"""
    print("\n=== Network Structure Optimization ===")
    
    # Calculate centrality measures
    try:
        # Degree centrality
        degree_cent = nx.degree_centrality(pkn_graph)
        
        # Betweenness centrality
        between_cent = nx.betweenness_centrality(pkn_graph)
        
        # Closeness centrality (for weakly connected components)
        close_cent = {}
        for component in nx.weakly_connected_components(pkn_graph):
            subgraph = pkn_graph.subgraph(component)
            if len(component) > 1:
                close_cent.update(nx.closeness_centrality(subgraph))
        
        # Find most important nodes
        print("Top 5 nodes by degree centrality:")
        top_degree = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:5]
        for node, cent in top_degree:
            print(f"  {node}: {cent:.3f}")
        
        print("\nTop 5 nodes by betweenness centrality:")
        top_between = sorted(between_cent.items(), key=lambda x: x[1], reverse=True)[:5]
        for node, cent in top_between:
            print(f"  {node}: {cent:.3f}")
            
        # Identify bottlenecks and key regulators
        bottlenecks = [node for node, cent in between_cent.items() if cent > 0.1]
        print(f"\nPotential bottleneck nodes: {bottlenecks}")
        
        return {
            'degree_centrality': degree_cent,
            'betweenness_centrality': between_cent,
            'closeness_centrality': close_cent,
            'bottlenecks': bottlenecks
        }
        
    except Exception as e:
        print(f"Error in structure optimization: {e}")
        return {}

# Method 4: Simulate network dynamics without experimental data
def simulate_network_dynamics(pkn_graph, boolean_rules):
    """Simulate Boolean network dynamics"""
    print("\n=== Boolean Network Simulation ===")
    
    try:
        nodes = list(pkn_graph.nodes())
        n_nodes = len(nodes)
        
        # Create random initial states for demonstration
        n_simulations = 5
        max_steps = 10
        
        for sim in range(n_simulations):
            print(f"\nSimulation {sim + 1}:")
            
            # Random initial state
            state = {node: np.random.choice([0, 1]) for node in nodes}
            print(f"Initial state: {state}")
            
            # Simulate dynamics
            trajectory = [state.copy()]
            
            for step in range(max_steps):
                new_state = state.copy()
                
                # Update based on simple rules (can be enhanced)
                for node in nodes:
                    predecessors = list(pkn_graph.predecessors(node))
                    if predecessors:
                        # Simple OR logic
                        new_state[node] = int(any(state[pred] for pred in predecessors))
                    # Input nodes keep their state or can be set externally
                
                # Check for steady state
                if new_state == state:
                    print(f"  Reached steady state at step {step + 1}: {new_state}")
                    break
                
                state = new_state
                trajectory.append(state.copy())
            
            if len(trajectory) == max_steps + 1:
                print(f"  Final state after {max_steps} steps: {state}")
        
    except Exception as e:
        print(f"Error in simulation: {e}")

# Main execution
def main():
    """Main function to demonstrate knowledge graph analysis"""
    print("Knowledge Graph Analysis without Experimental Data")
    print("=" * 50)
    
    # Load the PKN (Pathway Knowledge Network)
    G = load_sif("data/apoptosis.sif")
    # print(G.nodes())
    # print(G.edges(data=True))
    
    # Method 1: CellNOpt topology analysis
    sif_data = analyze_with_cellnopt(G)

    # Method 2: CASPO logic analysis
    boolean_rules = analyze_with_caspo(G)

    # Method 3: Network optimization
    centrality_metrics = optimize_network_structure(G)
    
    # Method 4: Simulate dynamics
    simulate_network_dynamics(G, boolean_rules)

    print("\n" + "=" * 50)
    print("Analysis complete! This demonstrates how to work with")
    print("knowledge graphs without experimental MIDAS data.")

if __name__ == "__main__":
    main()