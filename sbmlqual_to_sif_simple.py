#!/usr/bin/env python3
"""
SBMLqual to SIF Converter using existing packages

This script uses packages that are likely already available in your environment.
It provides multiple methods to convert SBMLqual to SIF format.

Requirements:
- xml.etree.ElementTree (built-in)
- networkx (optional)
"""

import xml.etree.ElementTree as ET
import argparse
from pathlib import Path
import sys

def extract_network_from_sbmlqual(sbml_file):
    """
    Extract network information from SBMLqual file using XML parsing.
    
    Returns:
        dict: Network information with nodes and edges
    """
    tree = ET.parse(sbml_file)
    root = tree.getroot()
    
    # Define namespaces - handle both prefixed and default namespaces
    qual_ns = 'http://www.sbml.org/sbml/level3/version1/qual/version1'
    core_ns = 'http://www.sbml.org/sbml/level3/version1/core'
    
    network = {
        'nodes': {},
        'edges': [],
        'metadata': {}
    }
    
    # Extract nodes (qualitative species) - try different approaches
    species_elements = []
    
    # Find all qualitative species elements
    for elem in root.iter():
        if elem.tag.endswith('qualitativeSpecies'):
            species_elements.append(elem)
    
    print(f"Debug: Found {len(species_elements)} qualitative species")
    
    for species in species_elements:
        # Handle namespaced attributes
        species_id = species.get(f'{{{qual_ns}}}id')
        species_name = species.get(f'{{{qual_ns}}}name', species_id)
        initial_level = species.get(f'{{{qual_ns}}}initialLevel', '0')
        max_level = species.get(f'{{{qual_ns}}}maxLevel', '1')
        
        if species_id:
            network['nodes'][species_id] = {
                'name': species_name,
                'id': species_id,
                'initialLevel': initial_level,
                'maxLevel': max_level
            }
            print(f"Debug: Added node {species_id}: {species_name}")
    
    # Extract edges (transitions)
    transition_elements = []
    
    # Find all transition elements
    for elem in root.iter():
        if elem.tag.endswith('transition'):
            transition_elements.append(elem)
    
    print(f"Debug: Found {len(transition_elements)} transitions")
    
    for transition in transition_elements:
        transition_id = transition.get(f'{{{qual_ns}}}id', f"transition_{len(network['edges'])}")
        
        # Find inputs in listOfInputs
        inputs = []
        for elem in transition.iter():
            if elem.tag.endswith('input'):
                source_id = elem.get(f'{{{qual_ns}}}qualitativeSpecies')
                sign = elem.get(f'{{{qual_ns}}}sign', 'unknown')
                threshold = elem.get(f'{{{qual_ns}}}thresholdLevel', '1')
                
                if source_id:
                    inputs.append({
                        'source': source_id,
                        'sign': sign,
                        'threshold': threshold
                    })
                    print(f"Debug: Found input {source_id} -> {sign}")
        
        # Find outputs in listOfOutputs
        outputs = []
        for elem in transition.iter():
            if elem.tag.endswith('output'):
                target_id = elem.get(f'{{{qual_ns}}}qualitativeSpecies')
                output_level = elem.get(f'{{{qual_ns}}}outputLevel', '1')
                
                if target_id:
                    outputs.append({
                        'target': target_id,
                        'level': output_level
                    })
                    print(f"Debug: Found output {target_id}")
        
        # Create edges
        for inp in inputs:
            for out in outputs:
                edge = {
                    'source': inp['source'],
                    'target': out['target'],
                    'sign': inp['sign'],
                    'threshold': inp['threshold'],
                    'output_level': out['level'],
                    'transition_id': transition_id
                }
                network['edges'].append(edge)
                print(f"Debug: Created edge {inp['source']} -> {out['target']} ({inp['sign']})")
    
    return network

def network_to_sif(network, interaction_types=None):
    """
    Convert network to SIF format.
    
    Args:
        network (dict): Network information
        interaction_types (dict): Mapping of signs to interaction types
    
    Returns:
        list: SIF formatted interactions
    """
    if interaction_types is None:
        interaction_types = {
            'positive': 'activates',
            'negative': 'inhibits',
            'unknown': 'regulates',
            'dual': 'regulates'
        }
    
    sif_lines = []
    
    for edge in network['edges']:
        source_id = edge['source']
        target_id = edge['target']
        sign = edge.get('sign', 'unknown')
        
        # Get node names
        source_name = network['nodes'].get(source_id, {}).get('name', source_id)
        target_name = network['nodes'].get(target_id, {}).get('name', target_id)
        
        # Determine interaction type
        interaction_type = interaction_types.get(sign, 'regulates')
        
        sif_line = f"{source_name}\t{interaction_type}\t{target_name}"
        sif_lines.append(sif_line)
    
    return sif_lines

def convert_sbmlqual_to_sif(input_file, output_file=None, verbose=True):
    """
    Main conversion function.
    """
    if verbose:
        print(f"Processing SBMLqual file: {input_file}")
    
    # Extract network
    try:
        network = extract_network_from_sbmlqual(input_file)
    except Exception as e:
        print(f"Error parsing SBMLqual file: {e}")
        return None
    
    if verbose:
        print(f"Found {len(network['nodes'])} nodes and {len(network['edges'])} edges")
    
    # Convert to SIF
    sif_interactions = network_to_sif(network)
    
    if verbose:
        print(f"Generated {len(sif_interactions)} SIF interactions")
    
    # Write output
    if output_file:
        with open(output_file, 'w') as f:
            f.write("\n".join(sif_interactions))
            f.write("\n")  # Add final newline
        if verbose:
            print(f"SIF file written to: {output_file}")
    
    return sif_interactions, network

def main():
    parser = argparse.ArgumentParser(description='Convert SBMLqual to SIF format')
    parser.add_argument('input', help='Input SBMLqual file (.xml)')
    parser.add_argument('output', nargs='?', help='Output SIF file (.sif)')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Verbose output')
    parser.add_argument('--show-network', action='store_true',
                       help='Show network information')
    
    args = parser.parse_args()
    
    # Determine output file
    if args.output is None:
        input_path = Path(args.input)
        args.output = input_path.with_suffix('.sif')
    
    # Convert
    try:
        result = convert_sbmlqual_to_sif(args.input, args.output, args.verbose)
        if result is None:
            sys.exit(1)
        
        sif_interactions, network = result
        
        if args.show_network:
            print("\\nNetwork Information:")
            print(f"Nodes ({len(network['nodes'])}):")
            for node_id, node_info in network['nodes'].items():
                print(f"  {node_id}: {node_info['name']}")
            
            print(f"\\nEdges ({len(network['edges'])}):")
            for edge in network['edges'][:10]:  # Show first 10
                print(f"  {edge['source']} -> {edge['target']} ({edge['sign']})")
            if len(network['edges']) > 10:
                print(f"  ... and {len(network['edges']) - 10} more")
        
        if args.verbose and sif_interactions:
            print("\\nFirst few SIF interactions:")
            for interaction in sif_interactions[:5]:
                print(f"  {interaction}")
            if len(sif_interactions) > 5:
                print(f"  ... and {len(sif_interactions) - 5} more")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
