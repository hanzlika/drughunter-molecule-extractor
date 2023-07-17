import argparse
import requests
from bs4 import BeautifulSoup
from calendar import month_name
from pdf2image import convert_from_bytes
from PIL import Image
import numpy as np
import os
from rdkit import Chem
import pandas as pd


def download_pdf(url:str):
    """
    Attempt to download PDF files from Drughunter Molecules of the Month based on the given URL.

    Params:
    - url (str): The URL of the web page to search for PDF files.

    Returns:
    - A list of tuples, where each tuple contains the file name and content of a PDF file found on the web page.
    - If no PDF files are found or an error occurs, an empty list is returned.
    """
    # headers simulating a browser request are necessary, without them the request status code returns 403 (forbidden)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # headers=headers prevents 403 forbidden response.status_code
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        # parses the web page for the first pdf file (assuming the molecules
        # of the month file is the first one on the web page)
        soup = BeautifulSoup(response.content, 'html.parser')
        pdf_links = soup.select('a[href$=".pdf"]')

        if pdf_links:
            pdf_files = []
            for pdf_link in pdf_links:
                file_name = pdf_link['href'].split('/')[-1]
                file_url = pdf_link['href']
            
                pdf_response = requests.get(file_url, headers=headers)

                if pdf_response.status_code == 200:
                    pdf_files.append((file_name, pdf_response.content))
                    print(file_name + " downloaded successfully.")
                else:
                    print("Failed to download the PDF:", file_name)
            return pdf_files
        else:
            print("PDF links not found on the page.")
    else:
        print("Failed to download the web page.")

    return []

def segment_pdf(pdf):
    from decimer_segmentation import segment_chemical_structures, segment_chemical_structures_from_file
    """
    Iterate through the PDF pages and use decimer_segmentation to get segmented images.

    Params:
    - pdf (bytes): The PDF file content in bytes.

    Returns:
    - List[Image]: A list of segmented images extracted from the PDF.
    """

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    segment_list = []
    pages = convert_from_bytes(pdf, 300)
    
    for page_number, page in enumerate(pages):

        # visualization can be set to True for a visual confirmation of the segmentation
        # expand=True yields better results than expand=False
        segments = segment_chemical_structures(np.array(page),
                                               expand=True,
                                               visualization=False)

        for segment in segments:
            image = Image.fromarray(segment)
            segment_list.append(image)

    return segment_list

def recognize_segments(image_list):
    from DECIMER import predict_SMILES
    """
    Iterate through the image list and attempt to recognize each image using DECIMER as a SMILES string.

    Params:
    - image_list (List[Image]): A list of images to recognize.

    Returns:
    - List[str]: A list of recognized SMILES strings with matching indexing.
    """
    print("Recognizing:")
    smiles_list = []
    for img in image_list:
        # Save the image as a temporary file
        img.save("tmp.png")

        # Use DECIMER to predict the SMILES string for the image
        # predict_SMILES function requires a path to the image as input
        SMILES = predict_SMILES("tmp.png")
        smiles_list.append(SMILES)

    # Remove the temporary image file
    if image_list:
        os.remove("tmp.png")

    return smiles_list

def smiles_to_inchi(smiles):
    """
    Convert a SMILES string to an InChI string representation.

    Params:
    - smiles (str): The SMILES string to convert.

    Returns:
    - str: The InChI string representation of the molecule.

    If the input SMILES string is not valid, it returns "Invalid SMILES".
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return "Invalid SMILES"

    inchi = Chem.MolToInchi(mol)
    return inchi

def smiles_list_to_inchi_list(smiles_list):
    """
    Convert a list of SMILES strings to a list of InChI string representations.

    Params:
    - smiles_list (List[str]): A list of SMILES strings to convert.

    Returns:
    - List[str]: A list of InChI string representations of the molecules.

    The function iterates through each SMILES string in the input list and uses the `smiles_to_inchi` function
    to convert it to the corresponding InChI string representation. The resulting InChI strings are appended
    to the `inchi_list` and returned at the end.
    """
    print("Converting smiles to inchi")
    inchi_list = []
    for smiles in smiles_list:
        inchi_list.append(smiles_to_inchi(smiles))
    return inchi_list
    
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

def export_results(segment_list, smiles_list, inchi_list, validated_list):
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
               "validated by Unichem": validated_list}

    # Create a DataFrame from the results dictionary
    df = pd.DataFrame(results)

    # Export the DataFrame to a CSV file
    df.to_csv('results.csv') 

def extract_molecules_from_url(url, can_choose):
    """
    Extract molecules from a given URL by downloading PDF files, segmenting them into images,
    recognizing SMILES strings, converting to InChI strings, and validating the InChI strings.

    Params:
    - url (str): The URL to download PDF files from.
    - can_choose (bool): Flag indicating if the user can choose a specific file to continue with.

    Returns:
    - Tuple[List[Tuple[str, Image]], List[str], List[str], List[bool]]: A tuple containing the segment list,
      SMILES list, InChI list, and validated list.

    The function attempts to download PDF files from the specified URL using the `download_pdf` function.
    It then prompts the user to choose a specific file if `can_choose` is `True`.
    After selecting the file(s), it segments the PDFs into images using the `segment_pdf` function.
    Then it recognizes the SMILES strings from the segmented images using the `recognize_segments` function.
    The recognized SMILES strings are converted to InChI strings using the `smiles_list_to_inchi_list` function.
    Finally, the InChI strings are validated using the `validate_inchi` function.

    The function returns the segment list, SMILES list, InChI list, and validated list as a tuple.
    """
    pdf_list = []
    
    print(f"Attempting to download pdf files from {url}")

    pdf_list += download_pdf(url)
    print(f"Downloaded {len(pdf_list)} PDF files from {url}.")

    for index, (filename, _) in enumerate(pdf_list):
        print(f"{index}: {filename}")

    # Prompt user to choose a filename
    if can_choose:
        while True:
            choice = input("Enter the index of the filename you want to continue with (or 'q' to quit): ")
            if choice.lower() == 'q':
                print("Quitting...")
                break
            try:
                index = int(choice)
                if -1 < index < len(pdf_list):
                    selected_files = [pdf_list[index]]
                    # Continue with the selected file
                    print("Selected file:", selected_files[0][0])
                    break
                else:
                    print("Invalid choice. Please try again.")
            except ValueError:
                print("Invalid input. Please enter a number or 'q' to quit.")
    else:
        selected_files = pdf_list

    segment_list = []
    for filename, pdf in selected_files:
        print(f"Attempting to segment {filename}")
        result = segment_pdf(pdf)
        to_append = [(filename, image) for image in result]
        print(f"Found {len(to_append)} segments in {filename}.")
        segment_list += to_append
    
    smiles_list = recognize_segments([segment for filename, segment in segment_list])

    inchi_list = smiles_list_to_inchi_list(smiles_list)

    validated_list = validate_inchi(inchi_list)

    return segment_list, smiles_list, inchi_list, validated_list

def extract_molecules_of_the_month(target_year):
    """
    Extract molecules of the month for a specific target year.

    Params:
    - target_year (int): The target year to extract molecules of the month.

    This function iterates over the months (from index 1 to 13) and constructs the URL for each month
    using the `target_year`. It then calls the `extract_molecules_from_url` function to extract molecules
    from each URL. The extracted segment, SMILES, InChI, and validation results are appended to their respective lists.

    Finally, the function calls the `export_results` function to export the results to a CSV file.
    """
    segment_list, smiles_list, inchi_list, validated_list = [], [], [], []
    for index in range(1, 13):
        url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{month_name[index].lower()}-{target_year}"
        tmp_segment_list, tmp_smiles_list, tmp_inchi_list, tmp_validated_list = extract_molecules_from_url(url, False)
        segment_list += tmp_segment_list
        smiles_list += tmp_smiles_list
        inchi_list += tmp_inchi_list
        validated_list += tmp_validated_list
    
    export_results(segment_list, smiles_list, inchi_list, validated_list)

def extract_molecules(url):
    """
    Extract molecules from a given URL.

    Params:
    - url (str): The URL to extract molecules from.

    This function calls the `extract_molecules_from_url` function to extract the segment list, SMILES list,
    InChI list, and validation list from the provided URL, allowing the user to choose a specific file.

    After extracting the molecules, the function calls the `export_results` function to export the results to a CSV file.
    """
    segment_list, smiles_list, inchi_list, validated_list = extract_molecules_from_url(url, True)

    export_results(segment_list, smiles_list, inchi_list, validated_list)


def main():
    """
    Main function to handle command-line arguments and execute the appropriate extraction based on the arguments.

    The function uses the `argparse` module to parse the command-line arguments.
    - If the `url` argument is provided, it calls the `extract_molecules` function with the provided URL.
    - If the `set` argument is not provided, it prints an error message.
    - If the `year` argument is not provided or is in an invalid format, it prints an error message.
    - If the `set` is "month", it calls the `extract_molecules_of_the_month` function with the provided year.
    - If the `set` is "drug", it calls the `extract_approved_drugs` function with the provided year.
    - If the `set` is not "month" or "drug", it prints an error message.

    Note: The code for the `extract_approved_drugs` function is not provided in the current context.
    """
    parser = argparse.ArgumentParser(description='DrugHunter extractor')
    parser.add_argument('-y', '--year', type=int, help='targeted year of drughunter sets')
    parser.add_argument('-s', '--set', type=str, help='targeted set (currently supports "month" for drugs of the month and \"drug\" for approval reviews from that year)')
    parser.add_argument('-u', '--url', type=str, help='url of webpage with targeted set (in case the format of drughunter url changes, which is likely)')
    args = parser.parse_args()
    
    # user requested a specific url to be checked, year and set is ignored
    if args.url:
        extract_molecules(args.url)
        return

    # user has not picked the set
    if not args.set:
        print("No set or url provided")
        return
    
    # user has not picked the year
    if not args.year:
        print("No year provided")
        return

    if len(str(args.year)) != 4:
        print("Invalid year format. Please provide a 4-digit year (YYYY).")
        return
    
    if args.set == "month":
        if (args.year < 2020):
            print("DrugHunter molecules of the month start at the year 2020, please provide a year equal or greater.")
            return
        extract_molecules_of_the_month(args.year)
        return

    if args.set == "drug":
        if (args.year < 2019):
            print("DrugHunter Drug approvals start at the year 2019, please provide a year equal or greater.")
        extract_approved_drugs(args.year)
        return

    print("The set parameter currently supports \"month\" and \"drug\" options, please use one of them, or provide a url to the web page you want pdf extracted from")

if __name__ == "__main__":
    main()