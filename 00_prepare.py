import pandas as pd
import numpy as np
import pyboolnet.file_exchange as FileExchange
import pyboolnet.interaction_graphs as IG
import pyboolnet.attractors as Attractors
import pyboolnet.state_space as StateSpace
import pyboolnet.state_transition_graphs as STGs
from pyboolnet.repository import get_primes
from pyboolnet.file_exchange import bnet2primes
from pyboolnet.prime_implicants import create_variables
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import re


# 1. Load your prior knowledge network
# Import BoolNet
boolnet = importr('BoolNet')

# Load SBML-qual file
r('net <- loadSBML("data/apoptosis.xml")')

# Print network info
r('print(net)')

# List all genes/nodes
genes = list(r('net$genes'))
print("Genes:", genes)

# Show transition functions (Boolean rules)
functions = r('net$functions')
print("Transition functions:")
for k, v in zip(functions.names, list(functions)):
    print(f"{k}: {v}")

# Find input nodes (cues/stimuli)
inputs = list(r('setdiff(net$genes, unlist(lapply(net$functions, all.vars)))'))
print("Inputs (cues/stimuli):", inputs)

# Find output nodes (readouts)
outputs = list(r('setdiff(net$genes, names(net$functions))'))
print("Outputs (readouts):", outputs)

# Find inhibitors: nodes that appear with a NOT (!) in any Boolean rule
inhibitors = set()
for rule in list(functions):
    found = re.findall(r'!([A-Za-z0-9_]+)', rule)
    inhibitors.update(found)
print("Inhibitors (appear as !X in rules):", sorted(inhibitors))

# Save as .bnet
r('saveNetwork(net, "data/apoptosis.bnet")')
print("Network saved as data/apoptosis.bnet")
primes = FileExchange.bnet2primes("data/apoptosis.bnet")

# 2. Generate synthetic data
# Boolean Network Rules + Perturbations → Simulated Dynamics → Synthetic Measurements
# Input: Gene regulatory rules from apoptosis.xml
# ↓
# Apply: Random perturbations (knockdowns, overexpressions)
# ↓
# Simulate: Network dynamics over time
# ↓
# Output: Time-series measurements (like MIDAS format)
def generate_synthetic_perturbations(primes, num_experiments=50):
    """Generate random perturbation experiments"""
    nodes = list(primes.keys())
    perturbations = []
    
    for i in range(num_experiments):
        # Random subset of nodes to perturb
        num_perturbed = np.random.randint(1, min(4, len(nodes)))
        perturbed_nodes = np.random.choice(nodes, num_perturbed, replace=False)
        
        perturbation = {}
        for node in nodes:
            if node in perturbed_nodes:
                perturbation[node] = np.random.choice([0, 1])
            else:
                perturbation[node] = None  # No perturbation
        
        perturbations.append(perturbation)
    
    return perturbations

def simulate_time_series(primes, perturbation, time_points=[0, 10, 30, 60]):
    """Simulate Boolean network dynamics under perturbation"""
    # Set initial state with perturbations
    initial_state = {}
    for node in primes.keys():
        if perturbation[node] is not None:
            initial_state[node] = perturbation[node]
        else:
            initial_state[node] = np.random.choice([0, 1])
    
    # Simulate dynamics
    states = [initial_state.copy()]
    current_state = initial_state.copy()
    
    for t in range(max(time_points)):
        # Update state using Boolean rules
        next_state = STGs.successor_synchronous(primes, current_state)
        states.append(next_state)
        current_state = next_state
    
    # Extract states at specified time points
    time_series = {}
    for tp in time_points:
        time_series[f'TR{tp}'] = states[min(tp, len(states)-1)]
    
    return time_series

def create_midas_format_data(primes, perturbations, time_points=[0, 10, 30, 60]):
    """Create synthetic data in MIDAS format"""
    all_data = []
    
    for exp_id, perturbation in enumerate(perturbations):
        time_series = simulate_time_series(primes, perturbation, time_points)
        
        for time_point in time_points:
            row = {'experiment': exp_id, 'time': time_point}
            
            # Add perturbation columns
            for node in primes.keys():
                if perturbation[node] is not None:
                    row[f'{node}i'] = perturbation[node]
                else:
                    row[f'{node}i'] = 0  # No perturbation
            
            # Add measurement columns
            state = time_series[f'TR{time_point}']
            for node in primes.keys():
                row[node] = state[node]
            
            all_data.append(row)
    
    return pd.DataFrame(all_data)

def add_noise_to_data(df, noise_level=0.1):
    """Add measurement noise to synthetic data"""
    measurement_cols = [col for col in df.columns if not col.endswith('i') and col not in ['experiment', 'time']]
    
    for col in measurement_cols:
        # Add Gaussian noise and clip to [0, 1]
        noise = np.random.normal(0, noise_level, len(df))
        df[col] = np.clip(df[col] + noise, 0, 1)
    
    return df


# Usage example:
if __name__ == "__main__":
    # Load your network
    primes = FileExchange.bnet2primes("data/apoptosis.bnet")

    # Generate perturbation experiments
    perturbations = generate_synthetic_perturbations(primes, num_experiments=50)
    
    # Create synthetic MIDAS-format data
    synthetic_data = create_midas_format_data(primes, perturbations)
    
    # Add measurement noise
    synthetic_data = add_noise_to_data(synthetic_data, noise_level=0.1)
    
    # Save to CSV
    synthetic_data.to_csv("data/synthetic_midas_data.csv", index=False)
    
    print(f"Generated synthetic data with {len(synthetic_data)} observations")
    print(f"Columns: {list(synthetic_data.columns)}")
    print(synthetic_data.head())