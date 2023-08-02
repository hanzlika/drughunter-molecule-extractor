from argparse import ArgumentParser
from numpy import asarray
from calendar import month_name
from pdf_extraction.pdf_extraction import download_pdf
from segmentation.segment_pdf import segment_pdf
from recognition.recognize_segments_molscribe import recognize_segments
from conversion.tools import smiles_list_to_inchikey_list
from validation.validate_with_chembl_webresource import validate_inchikey_list
from export.export_results import export_to_csv

def extract_molecules_from_pdfs(pdfs):
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
    extraction_results = {}


    # segmentation_results -> list of tuples (source filename, segment)
    segmentation_results = segment_pdf(pdfs, save_to_directory=True, directory="segmentation/segments") 

    segments = [segment for (_, segment) in segmentation_results]
    extraction_results['source'] = [source for source, _ in segmentation_results]
    
    extraction_results.update(recognize_segments(segments))

    extraction_results['validation'] = validate_inchikey_list(extraction_results.get('inchikey'))

    if extraction_results:
        export_to_csv(extraction_results)


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

    pdfs = []
    for index in range(1, 13):
        url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{month_name[index].lower()}-{target_year}"
        pdfs += download_pdf(url, download_all=True)

    extract_molecules_from_pdfs(pdfs)


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
    parser = ArgumentParser(description='DrugHunter extractor')
    parser.add_argument('-y', '--year', type=int, help='targeted year of drughunter sets')
    parser.add_argument('-s', '--set', type=str, help='targeted set (currently supports "month" for drugs of the month and \"drug\" for approval reviews from that year)')
    parser.add_argument('-u', '--url', type=str, help='url of webpage with targeted set (in case the format of drughunter url changes, which is likely)')
    args = parser.parse_args()
    
    # user requested a specific url to be checked, year and set is ignored
    if args.url:
        extract_molecules_from_pdfs(download_pdf(args.url, download_all=False))
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