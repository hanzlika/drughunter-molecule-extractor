import pandas as pd

def export_results(segment_list, smiles_list, inchi_list, inchikey_list, validated_list, validated_inchi_key_list):
    """
    Export the segmentation, SMILES, InChI, and validation results to a CSV file.

    Params:
    - segment_list (List[Tuple[str, bytes]]): A list of tuples containing the filename and segmented image data.
    - smiles_list (List[str]): A list of SMILES strings.
    - inchi_list (List[str]): A list of InChI strings.
    - validated_list (List[bool]): A list of boolean values indicating the validation status of each InChI.

    The function creates a dictionary named `results` with four keys: "source", "smiles", "inchi", and "validated by Unichem".
    The "source" key corresponds to the filenames extracted from the `segment_list`.
    The "smiles" key corresponds to the `smiles_list`.
    The "inchi" key corresponds to the `inchi_list`.
    The "validated by Unichem" key corresponds to the `validated_list`.
    
    The `results` dictionary is then used to create a DataFrame using pandas, and the DataFrame is exported to a CSV file named 'results.csv'.
    """
    print("Exporting results.")    

    # Create a dictionary with the results
    results = {"source": [filename for filename, _ in segment_list],
               "smiles": smiles_list, 
               "inchi": inchi_list,
               "inchikey": inchikey_list, 
               "inchi found in unichem": validated_list,
               "inchikey found in unichem connectivity": validated_inchi_key_list}

    # Create a DataFrame from the results dictionary
    df = pd.DataFrame(results)

    # Export the DataFrame to a CSV file
    df.to_csv('results_with_validation_comparison.csv') 