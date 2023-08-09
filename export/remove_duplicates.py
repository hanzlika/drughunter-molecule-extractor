def remove_duplicates(dict_with_duplicates : dict) -> dict:
    """
    Removes rows of information from a dict so that no duplicate inchi keys are present

    Params:
        - dict_with_duplicates (dict): dict unfiltered for duplicate information

    Returns:
        - dict_without_duplicates (dict): dict filtered for duplicate information
    """
    inchikeys = set()
    dict_without_duplicates = {key: [] for key in dict_with_duplicates.keys()}
    removed_count = 0
    for index, inchikey in enumerate(dict_with_duplicates['inchikey']):
        if inchikey not in inchikeys:
            for key in dict_without_duplicates.keys():
                dict_without_duplicates.get(key).append((dict_with_duplicates.get(key))[index])
            if inchikey:
                inchikeys.add(inchikey)
        else:
            removed_count += 1
    print(f"{removed_count} rows were removed due to being duplicates of previous rows")
    return dict_without_duplicates