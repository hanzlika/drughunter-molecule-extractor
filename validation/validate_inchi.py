import requests

def validate_inchi(inchi_list):
    """
    Validate InChI strings through the UNICHEM API.

    Params:
    - inchi_list (List[str]): A list of InChI strings to validate.

    Returns:
    - List[bool]: A list of boolean values indicating the validation status for each InChI string.

    The function iterates through each InChI string in the input list and sends a POST request to the UNICHEM API
    to validate the InChI. The response from the API is checked, and if the response is not "Not found",
    indicating that the InChI is valid, the function appends `True` to the `validated_list`.
    Otherwise, it appends `False` to the `validated_list`.

    The resulting `validated_list` is returned at the end.
    """
    print("Validating inchi through UNICHEM.")

    # update in case the api changes again
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

    payload = {"type": "inchi", "compound": ""}

    headers = {"Content-Type": "application/json"}

    validated_list = []
    for inchi in inchi_list:
        payload["compound"] = inchi
        response = requests.request("POST", url, json=payload, headers=headers)
        resp_dict = response.json()
        if resp_dict["response"] != "Not found":
            validated_list.append(True)
        else:
            validated_list.append(False)

    return validated_list