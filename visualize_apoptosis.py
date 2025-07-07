#!/usr/bin/env python3
"""
Alternative visualization script for apoptosis.xml using Python libraries
This script extracts network information and creates a visualization using networkx and matplotlib
"""

import os
import xml.etree.ElementTree as ET
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from cellnopt2py import CellNOptAnalyzer

def parse_sbml_to_networkx(xml_file):
    """
    Parse SBML file and create a NetworkX graph
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Find the namespace - try different common namespaces
        ns_variations = [
            {'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
             'qual': 'http://www.sbml.org/sbml/level3/version1/qual/version1'},
            {'sbml': 'http://www.sbml.org/sbml/level2/version4',
             'qual': 'http://www.sbml.org/sbml/level3/version1/qual/version1'},
            {}  # No namespace
        ]
        
        G = nx.DiGraph()
        
        for ns in ns_variations:
            try:
                # Find qualitative species (nodes)
                if ns:
                    species_elements = root.findall('.//qual:qualitativeSpecies', ns)
                    if not species_elements:
                        species_elements = root.findall('.//qualitativeSpecies')
                else:
                    species_elements = root.findall('.//qualitativeSpecies')
                
                for species in species_elements:
                    species_id = species.get('id')
                    species_name = species.get('name')
                    if species_id:  # Only add if id exists
                        node_label = species_name if species_name else species_id
                        G.add_node(species_id, name=node_label, label=node_label)
                
                # Find transitions (edges/reactions)
                if ns:
                    transitions = root.findall('.//qual:transition', ns)
                    if not transitions:
                        transitions = root.findall('.//transition')
                else:
                    transitions = root.findall('.//transition')
                
                for transition in transitions:
                    # Get inputs and outputs
                    if ns:
                        inputs = transition.findall('.//qual:input', ns)
                        outputs = transition.findall('.//qual:output', ns)
                        if not inputs:
                            inputs = transition.findall('.//input')
                        if not outputs:
                            outputs = transition.findall('.//output')
                    else:
                        inputs = transition.findall('.//input')
                        outputs = transition.findall('.//output')
                    
                    for input_elem in inputs:
                        input_species = input_elem.get('qualitativeSpecies')
                        for output_elem in outputs:
                            output_species = output_elem.get('qualitativeSpecies')
                            if input_species and output_species and input_species in G.nodes() and output_species in G.nodes():
                                G.add_edge(input_species, output_species)
                
                if len(G.nodes()) > 0:  # If we found nodes, break
                    break
                    
            except Exception as e:
                print(f"Trying namespace variant: {e}")
                continue
        
        # If still no nodes found, try to parse any species elements
        if len(G.nodes()) == 0:
            print("Trying fallback parsing...")
            all_elements = root.findall('.//*')
            for elem in all_elements:
                if 'species' in elem.tag.lower() or 'node' in elem.tag.lower():
                    elem_id = elem.get('id')
                    elem_name = elem.get('name')
                    if elem_id:
                        node_label = elem_name if elem_name else elem_id
                        G.add_node(elem_id, name=node_label, label=node_label)
        
        return G if len(G.nodes()) > 0 else None
    
    except Exception as e:
        print(f"Error parsing SBML file: {e}")
        return None

def plot_network_with_python(G, output_file="output/apoptosis_network_python.png"):
    """
    Create a network plot using matplotlib and networkx
    """
    if len(G.nodes()) == 0:
        print("No nodes to plot!")
        return
    
    plt.figure(figsize=(16, 12))
    
    # Use hierarchical layout for better visualization of biological networks
    try:
        # Try to use graphviz layout if available
        pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
    except:
        try:
            # Fallback to spring layout with more iterations
            pos = nx.spring_layout(G, k=3, iterations=100, seed=42)
        except:
            # Last resort: random layout
            pos = nx.random_layout(G, seed=42)
    
    # Create node colors based on in/out degree
    node_colors = []
    for node in G.nodes():
        in_degree = G.in_degree(node)
        out_degree = G.out_degree(node)
        if in_degree == 0:
            node_colors.append('lightgreen')  # Source nodes
        elif out_degree == 0:
            node_colors.append('lightcoral')  # Sink nodes
        else:
            node_colors.append('lightblue')   # Intermediate nodes
    
    # Get node labels
    labels = {}
    for node in G.nodes():
        if 'label' in G.nodes[node]:
            labels[node] = G.nodes[node]['label']
        elif 'name' in G.nodes[node]:
            labels[node] = G.nodes[node]['name']
        else:
            labels[node] = str(node)
    
    # Draw the network
    nx.draw(G, pos, 
            labels=labels,
            node_color=node_colors,
            node_size=2000,
            font_size=10,
            font_weight='bold',
            arrows=True,
            arrowsize=25,
            edge_color='gray',
            linewidths=2,
            alpha=0.9)
    
    # Add a legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgreen', 
                   markersize=10, label='Source nodes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightcoral', 
                   markersize=10, label='Sink nodes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue', 
                   markersize=10, label='Intermediate nodes')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.title(f"Apoptosis Network Model\n({len(G.nodes())} nodes, {len(G.edges())} edges)", 
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Network plot saved to {output_file}")
    
    # Also save network statistics
    stats_file = output_file.replace('.png', '_stats.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Network Statistics\n")
        f.write(f"==================\n\n")
        f.write(f"Number of nodes: {len(G.nodes())}\n")
        f.write(f"Number of edges: {len(G.edges())}\n")
        f.write(f"Average degree: {sum(dict(G.degree()).values()) / len(G.nodes()):.2f}\n\n")
        f.write(f"Nodes:\n")
        for node in sorted(G.nodes()):
            f.write(f"  {node} (in: {G.in_degree(node)}, out: {G.out_degree(node)})\n")
        f.write(f"\nEdges:\n")
        for edge in sorted(G.edges()):
            f.write(f"  {edge[0]} -> {edge[1]}\n")
    print(f"Network statistics saved to {stats_file}")

def extract_network_with_r(xml_file):
    """
    Use R to extract network information and return as Python data
    """
    try:
        analyzer = CellNOptAnalyzer()
        analyzer.load_network(xml_file, USER_SIF=False)
        
        # Extract network information using R
        import rpy2.robjects as robjects
        import numpy as np
        
        # Get species names
        species = robjects.r('pknmodel$namesSpecies')
        species_list = [str(s) for s in species]
        
        # Get reaction matrix - convert to numpy array for easier handling
        reac_mat_r = robjects.r('pknmodel$interMat')
        reac_mat = np.array(reac_mat_r)
        
        # Convert to NetworkX graph
        G = nx.DiGraph()
        
        # Add nodes
        for species_name in species_list:
            G.add_node(species_name)
        
        # Add edges based on interaction matrix
        if reac_mat.ndim == 2:
            for i in range(reac_mat.shape[0]):  # For each reaction
                inputs = []
                outputs = []
                for j in range(reac_mat.shape[1]):  # For each species
                    if reac_mat[i, j] == -1:  # Input
                        inputs.append(species_list[j])
                    elif reac_mat[i, j] == 1:  # Output
                        outputs.append(species_list[j])
                
                # Create edges from inputs to outputs
                for input_node in inputs:
                    for output_node in outputs:
                        G.add_edge(input_node, output_node)
        
        return G
    
    except Exception as e:
        print(f"Error extracting network with R: {e}")
        return None

def main():
    """Main function to visualize apoptosis network"""
    
    xml_file = "data/apoptosis.xml"
    
    if not os.path.exists(xml_file):
        print(f"Error: {xml_file} not found")
        return
    
    os.makedirs("output", exist_ok=True)
    
    print("Attempting to visualize apoptosis network...")
    
    # Method 1: Try R-based extraction
    print("\nMethod 1: Using R to extract network...")
    G_r = extract_network_with_r(xml_file)
    if G_r and len(G_r.nodes()) > 0:
        plot_network_with_python(G_r, "output/apoptosis_network_from_r.png")
        print(f"Network has {len(G_r.nodes())} nodes and {len(G_r.edges())} edges")
    
    # Method 2: Try direct XML parsing
    print("\nMethod 2: Direct XML parsing...")
    G_xml = parse_sbml_to_networkx(xml_file)
    if G_xml and len(G_xml.nodes()) > 0:
        plot_network_with_python(G_xml, "output/apoptosis_network_from_xml.png")
        print(f"Network has {len(G_xml.nodes())} nodes and {len(G_xml.edges())} edges")
    
    print("\nDone! Check the output/ directory for network visualizations.")

if __name__ == "__main__":
    main()
