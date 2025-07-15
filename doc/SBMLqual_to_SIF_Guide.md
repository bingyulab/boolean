# SBMLqual to SIF Conversion Guide

## Overview

Simple Interaction Format (SIF) is a tab-delimited network format where each line represents an interaction:
```
source_node    interaction_type    target_node
```

## Usage in Your Boolean Modeling Pipeline
- **`sbmlqual_to_sif_simple.py`** - Pure Python solution using built-in XML parser
- **`sbmlqual_to_sif.py`** - Advanced solution using libSBML (requires installation)


1. **Convert SBMLqual to SIF:**
   ```bash
   python sbmlqual_to_sif_simple.py model.xml model.sif
   ```

2. **Use SIF in analysis tools:**
   - Import into Cytoscape for visualization
   - Convert to other Boolean network formats
   - Use with PyBoolNet for analysis
   - Input for network analysis in R/Python

3. **Integration with your workflow:**
   - SIF format is widely supported
   - Can be converted to adjacency matrices
   - Compatible with most network analysis tools

## Advanced Options

### Custom Interaction Types
```bash
python sbmlqual_to_sif_simple.py model.xml model.sif --interaction-type "regulates"
```

### Verbose Output
```bash
python sbmlqual_to_sif_simple.py model.xml --verbose --show-network
```

### Network Information
The converter also extracts additional network metadata:
- Node properties (maxLevel, initialLevel)
- Transition details (thresholds, function terms)
- Can be extended for more complex conversions

## Recommendations

1. **For simple conversions:** Use `sbmlqual_to_sif_simple.py` (no dependencies)
2. **For complex models:** Install libSBML and use `sbmlqual_to_sif.py`
3. **For interactive analysis:** Use GINsim GUI tool
4. **For R workflows:** Use BoolNet package
5. **For Python workflows:** Combine with PyBoolNet or NetworkX
