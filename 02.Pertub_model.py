import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.packages import importr


def get_descendants(mat, start_idx, visited=None):
    if visited is None:
        visited = set()
    for j, val in enumerate(mat[start_idx]):
        if val != 0 and j not in visited and j < mat.shape[0]:
            visited.add(j)
            get_descendants(mat, j, visited)
    return visited

def modify_model(change_percent):
    """
    Modify the model by perturbing a percentage of nodes.
    Args:
        change_percent: Percentage of nodes to perturb (0.0 to 1.0).
    Returns:
        None, modifies the model in place.
    """
    # Ensure change_percent is between 0 and 1
    if not (0 <= change_percent <= 1):
        raise ValueError("change_percent must be between 0 and 1")
    model_size = int(r('nrow(model$interMat)')[0])
    num_perturb = max(1, int(model_size * change_percent))

    # Get node names
    node_names = list(r('rownames(model$interMat)'))

    # Randomly select nodes to perturb
    perturb_nodes = np.random.choice(node_names, num_perturb, replace=False)

    for node in perturb_nodes:
        x = np.random.uniform(0, 3)
        print(f"Perturbing node: {node}, random x={x:.2f}")
        if x < 1:
            # Operator insertion: randomly connect to another node with AND/OR/NOT
            target = np.random.choice([n for n in node_names if n != node])
            op = np.random.choice(['AND', 'OR', 'NOT'])
            print(f"  Operator insertion: {op} between {node} and {target}")
            idx_node = node_names.index(node) + 1
            idx_target = node_names.index(target) + 1
            if op == 'NOT':
                r(f"model$notMat[{idx_node},{idx_target}] <- 1")
            else:
                if op == 'AND':
                    r(f"model$interMat[{idx_node},{idx_target}] <- 1")
                    r(f"model$interMat[{idx_target},{idx_node}] <- -1")
                else:
                    r(f"model$interMat[{idx_node},{idx_target}] <- -1")
                    r(f"model$interMat[{idx_target},{idx_node}] <- 1")
        elif 1 <= x and x < 2:
            # Subtree deletion: randomly zero out a row/column
            print(f"  Subtree deletion at {node}")
            idx_node = node_names.index(node)
            # Get interMat as numpy array
            intermat = np.array(r('as.matrix(model$interMat)'))
            # Find descendants using DFS
            descendants = get_descendants(intermat, idx_node)
            # Include the node itself
            to_delete = set([idx_node]) | descendants
            for idx in to_delete:
                r(f"model$interMat[{idx+1},] <- 0")
                r(f"model$interMat[,{idx+1}] <- 0")
                r(f"model$notMat[{idx+1},] <- 0")
                r(f"model$notMat[,{idx+1}] <- 0")
            print(f"    Deleted subtree nodes: {[node_names[i] for i in to_delete]}")

        else:
            # Operator exchange: swap two nonzero entries in interMat or notMat
            print(f"  Operator exchange")
            intermat = np.array(r('as.matrix(model$interMat)'))
            notmat = np.array(r('as.matrix(model$notMat)'))
            # Find nonzero indices
            nz_inter = np.argwhere(intermat != 0)
            nz_not = np.argwhere(notmat != 0)
            if len(nz_inter) > 1:
                idx1, idx2 = np.random.choice(len(nz_inter), 2, replace=False)
                i1, j1 = nz_inter[idx1]
                i2, j2 = nz_inter[idx2]
                print(f"    Swapping interMat[{i1+1},{j1+1}] <-> interMat[{i2+1},{j2+1}]")
                val1 = intermat[i1, j1]
                val2 = intermat[i2, j2]
                r(f"tmp <- model$interMat[{i1+1},{j1+1}]")
                r(f"model$interMat[{i1+1},{j1+1}] <- model$interMat[{i2+1},{j2+1}]")
                r(f"model$interMat[{i2+1},{j2+1}] <- tmp")
            elif len(nz_not) > 1:
                idx1, idx2 = np.random.choice(len(nz_not), 2, replace=False)
                i1, j1 = nz_not[idx1]
                i2, j2 = nz_not[idx2]
                print(f"    Swapping notMat[{i1+1},{j1+1}] <-> notMat[{i2+1},{j2+1}]")
                r(f"tmp <- model$notMat[{i1+1},{j1+1}]")
                r(f"model$notMat[{i1+1},{j1+1}] <- model$notMat[{i2+1},{j2+1}]")
                r(f"model$notMat[{i2+1},{j2+1}] <- tmp")
            else:
                print("    No operator to exchange.")

    # The modified model is now in R as 'model'
    # You can now use it for further CellNOptR/MEIGO/Caspo analysis
    # Save the modified model as an RData file
    r(f'save(model, file="output/ModifiedToyModel_{change_percent:.1%}.RData")')
    r(f'writeSIF(model, "output/caspo/ModifiedToyModel_{change_percent:.1%}.sif")')
    print("Model perturbation complete.")
    
    
if __name__ == "__main__":
    # Load the CellNOptR package
    cellnopt = importr('CellNOptR')

    # Load the ToyModel data
    r('data("ToyModel", package="CellNOptR")')
    r('model <- ToyModel')

    # Parameters
    for change_percent in np.linspace(0, 1, 11)[1:-1]:
        print(f"Modifying model with {change_percent:.1%} perturbation...")
        # Modify the model with a x% perturbation
        modify_model(change_percent=change_percent)

    # Load the modified model back into Python if needed
    # r(f'load("output/ModifiedToyModel_{change_percent:.1%}.RData")')