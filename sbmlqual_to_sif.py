#!/usr/bin/env python3
"""
SBMLqual to SIF Converter

This script converts SBML qualitative models (SBMLqual) to Simple Interaction Format (SIF).
SIF format: source_node <interaction_type> target_node

Requirements:
- python-libsbml
- networkx (optional, for network analysis)

Usage:
    python sbmlqual_to_sif.py input.xml output.sif
"""

import sys
import argparse
from pathlib import Path

try:
    import libsbml
    LIBSBML_AVAILABLE = True
except ImportError:
    LIBSBML_AVAILABLE = False
    print("Warning: libsbml not available. Install with: pip install python-libsbml")

def parse_sbmlqual_to_sif(sbml_file, output_file=None, interaction_type="regulates"):
    """
    Convert SBMLqual file to SIF format.
    
    Args:
        sbml_file (str): Path to SBMLqual file
        output_file (str): Path to output SIF file (optional)
        interaction_type (str): Default interaction type for SIF
    
    Returns:
        list: List of SIF interactions
    """
    if not LIBSBML_AVAILABLE:
        raise ImportError("libsbml is required. Install with: pip install python-libsbml")
    
    # Read SBML document
    document = libsbml.readSBMLFromFile(str(sbml_file))
    
    if document.getNumErrors() > 0:
        print(f"Errors reading SBML file {sbml_file}:")
        document.printErrors()
        return []
    
    model = document.getModel()
    if model is None:
        print("No model found in SBML file")
        return []
    
    # Check if qual plugin is available
    qual_plugin = model.getPlugin("qual")
    if qual_plugin is None:
        print("No qualitative models plugin found")
        return []
    
    sif_interactions = []
    
    # Extract qualitative species (nodes)
    species_dict = {}
    for i in range(qual_plugin.getNumQualitativeSpecies()):
        qs = qual_plugin.getQualitativeSpecies(i)
        species_dict[qs.getId()] = qs.getName() if qs.getName() else qs.getId()
    
    # Extract transitions (interactions)
    for i in range(qual_plugin.getNumTransitions()):
        transition = qual_plugin.getTransition(i)
        
        # Get outputs (target nodes)
        outputs = []
        for j in range(transition.getNumOutputs()):
            output = transition.getOutput(j)
            target_id = output.getQualitativeSpecies()
            target_name = species_dict.get(target_id, target_id)
            outputs.append(target_name)
        
        # Get inputs (source nodes)
        inputs = []
        for j in range(transition.getNumInputs()):
            input_element = transition.getInput(j)
            source_id = input_element.getQualitativeSpecies()
            source_name = species_dict.get(source_id, source_id)
            
            # Determine interaction type based on input properties
            input_type = interaction_type
            if input_element.isSetSign():
                sign = input_element.getSign()
                if sign == libsbml.INPUT_SIGN_POSITIVE:
                    input_type = "activates"
                elif sign == libsbml.INPUT_SIGN_NEGATIVE:
                    input_type = "inhibits"
                else:
                    input_type = "regulates"
            
            inputs.append((source_name, input_type))
        
        # Create SIF interactions
        for source, int_type in inputs:
            for target in outputs:
                sif_interactions.append(f"{source}\t{int_type}\t{target}")
    
    # Write to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write("\n".join(sif_interactions)); f.write("\n")
        print(f"SIF file written to: {output_file}")
    
    return sif_interactions

def parse_sbmlqual_simple(sbml_file, output_file=None):
    """
    Simple parser for SBMLqual without libsbml (XML parsing).
    This is a basic implementation that extracts network structure.
    """
    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        raise ImportError("XML parsing not available")
    
    # Parse XML
    tree = ET.parse(sbml_file)
    root = tree.getroot()
    
    # Find SBML namespace
    ns = {'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
          'qual': 'http://www.sbml.org/sbml/level3/version1/qual/version1'}
    
    sif_interactions = []
    
    # Find qualitative species
    species_dict = {}
    for qs in root.findall('.//qual:qualitativeSpecies', ns):
        species_id = qs.get('id')
        species_name = qs.get('name', species_id)
        species_dict[species_id] = species_name
    
    # Find transitions
    for transition in root.findall('.//qual:transition', ns):
        # Get inputs
        inputs = []
        for inp in transition.findall('.//qual:input', ns):
            source_id = inp.get('qualitativeSpecies')
            source_name = species_dict.get(source_id, source_id)
            sign = inp.get('sign', 'unknown')
            
            if sign == 'positive':
                int_type = 'activates'
            elif sign == 'negative':
                int_type = 'inhibits'
            else:
                int_type = 'regulates'
            
            inputs.append((source_name, int_type))
        
        # Get outputs
        outputs = []
        for out in transition.findall('.//qual:output', ns):
            target_id = out.get('qualitativeSpecies')
            target_name = species_dict.get(target_id, target_id)
            outputs.append(target_name)
        
        # Create SIF interactions
        for source, int_type in inputs:
            for target in outputs:
                sif_interactions.append(f"{source}\t{int_type}\t{target}")
    
    # Write to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write("\n".join(sif_interactions)); f.write("\n")
        print(f"SIF file written to: {output_file}")
    
    return sif_interactions

def main():
    parser = argparse.ArgumentParser(description='Convert SBMLqual to SIF format')
    parser.add_argument('input', help='Input SBMLqual file (.xml)')
    parser.add_argument('output', nargs='?', help='Output SIF file (.sif)')
    parser.add_argument('--simple', action='store_true', 
                       help='Use simple XML parser instead of libsbml')
    parser.add_argument('--interaction-type', default='regulates',
                       help='Default interaction type (default: regulates)')
    
    args = parser.parse_args()
    
    # Determine output file
    if args.output is None:
        input_path = Path(args.input)
        args.output = input_path.with_suffix('.sif')
    
    # Convert
    try:
        if args.simple or not LIBSBML_AVAILABLE:
            print("Using simple XML parser...")
            interactions = parse_sbmlqual_simple(args.input, args.output)
        else:
            print("Using libsbml parser...")
            interactions = parse_sbmlqual_to_sif(args.input, args.output, args.interaction_type)
        
        print(f"Extracted {len(interactions)} interactions")
        if len(interactions) > 0:
            print("First few interactions:")
            for i, interaction in enumerate(interactions[:5]):
                print(f"  {interaction}")
            if len(interactions) > 5:
                print(f"  ... and {len(interactions) - 5} more")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
