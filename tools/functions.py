import re
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def parse_namedlist_to_dataframes(named_list):
    """
    Parse rpy2 NamedList object containing attractor comparison results
    """
    results = {}
    
    # Get the names (keys) from the NamedList
    names = list(named_list.names)
    print(f"Found {len(names)} files: {names}")
    
    for name in names:
        print(f"\nProcessing: {name}")
        
        # Access the R DataFrame for this name
        r_dataframe = named_list[name]
        print(f"R DataFrame type: {type(r_dataframe)}")
        
        # Convert R DataFrame to pandas DataFrame
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.conversion.rpy2py(r_dataframe)
                results[name] = df
                print(f"Successfully converted {name} to pandas DataFrame")
                print(f"Shape: {df.shape}")
                print(f"Columns: {list(df.columns)}")
                print(f"Index: {list(df.index)}")
        except Exception as e:
            print(f"Conversion failed for {name}: {e}")
            # Try manual conversion
            df = manual_convert_r_dataframe(r_dataframe)
            results[name] = df
            print(f"Manual conversion successful for {name}")
    
    return results


def parse_rpy2_results(r_list_vector):
    """
    Parse RPy2 ListVector containing attractor comparison results using pandas2ri
    """
    results = {}
    
    # Enable pandas conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        # Convert the R list to Python dict
        named_list = ro.conversion.rpy2py(r_list_vector)
        names = list(named_list.names())
        print(f"Found {len(names)} files: {names}")
        
        for i, name in enumerate(names):
            print(f"\nProcessing: {name}")
            
            # Access the R DataFrame for this name
            r_dataframe = named_list[i]
            print(f"R DataFrame type: {type(r_dataframe)}")
            
            df = ro.conversion.rpy2py(r_dataframe)
            results[name] = df
            print(f"Successfully converted {name} to pandas DataFrame")
    return results

# Convert R list to Python dict
def rlist_to_pydict(rlist):
    py_dict = {}
    for name in rlist.names:
        value = rlist.rx2(name)
        # Convert R vectors to Python lists
        if hasattr(value, 'tolist'):
            py_dict[name] = value.tolist()
        else:
            py_dict[name] = value
    return py_dict
    
    
def analyze_attractor_performance(df):
    """
    Analyze attractor performance metrics
    """
    print("\nPerformance Analysis:")
    print(f"Best accuracy: {df['accuracy'].max():.4f} (Attractor: {df['accuracy'].idxmax()})")
    print(f"Best F1 score: {df['f1'].max():.4f} (Attractor: {df['f1'].idxmax()})")
    print(f"Average hamming distance: {df['hamming'].mean():.2f}")
    print(f"Average normalized hamming: {df['norm_hamming'].mean():.4f}")
    
    # Count valid precision/recall values
    valid_precision = df['precision'].notna().sum()
    valid_f1 = df['f1'].notna().sum()
    print(f"Attractors with valid precision: {valid_precision}/{len(df)}")
    print(f"Attractors with valid F1: {valid_f1}/{len(df)}")