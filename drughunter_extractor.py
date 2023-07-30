print("Importing (decimer may take up to a minute to import)")

import argparse
import requests
from calendar import month_name
from PIL import Image
import os
import pandas as pd
from pdf_extraction.extract_pdf import download_pdf
from segmentation.segment_pdf import segment_pdf
from recognition.recognize_segments import recognize_segments
from conversion.tools import smiles_list_to_inchi_list, smiles_list_to_inchikey_list
from validation.validate_inchi import validate_inchi
from validation.validate_with_chembl_webresource import validate_inchikey_list
from export.export_results import export_results

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
    if not pdf_list:
        return [], [], [], []

    for index, (filename, _) in enumerate(pdf_list):
        print(f"{index}: {filename}")

    # Prompt user to choose a filename
    if can_choose:
        while True:
            choice = input("Enter the index of the filename you want to continue with (or 'q' to quit): ")
            if choice.lower() == 'q':
                print("Quitting...")
                return [],[],[],[]
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
    inchikey_list = smiles_list_to_inchikey_list(smiles_list)

    validated_list = validate_inchi(inchi_list)
    validated_inchi_key_list = validate_inchikey_list(inchikey_list)

    return segment_list, smiles_list, inchi_list, inchikey_list, validated_list, validated_inchi_key_list

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
    segment_list, smiles_list, inchi_list, inchikey_list,  validated_list, validated_inchi_key_list = [], [], [], [], [], []
    for index in range(1, 13):
        url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{month_name[index].lower()}-{target_year}"
        tmp_segment_list, tmp_smiles_list, tmp_inchi_list, tmp_inchikey_list, tmp_validated_list, tmp_validated_inchi_key_list = extract_molecules_from_url(url, False)
        segment_list += tmp_segment_list
        smiles_list += tmp_smiles_list
        inchi_list += tmp_inchi_list
        inchikey_list += tmp_inchikey_list
        validated_list += tmp_validated_list
        validated_inchi_key_list += tmp_validated_inchi_key_list
    
    if segment_list and smiles_list and inchi_list and validated_list:
        export_results(segment_list, smiles_list, inchi_list, inchikey_list, validated_list, validated_inchi_key_list)

def extract_molecules(url):
    """
    Extract molecules from a given URL.

    Params:
    - url (str): The URL to extract molecules from.

    This function calls the `extract_molecules_from_url` function to extract the segment list, SMILES list,
    InChI list, and validation list from the provided URL, allowing the user to choose a specific file.

    After extracting the molecules, the function calls the `export_results` function to export the results to a CSV file.
    """
    segment_list, smiles_list, inchi_list, inchikey_list, validated_list, validated_inchi_key_list = extract_molecules_from_url(url, True)

    if segment_list and smiles_list and inchi_list and validated_list:
        export_results(segment_list, smiles_list, inchi_list, inchikey_list, validated_list, validated_inchi_key_list)


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