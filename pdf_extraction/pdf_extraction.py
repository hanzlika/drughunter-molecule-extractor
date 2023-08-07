import requests
import os
from bs4 import BeautifulSoup


def filter_pdf_links_by_choice(pdf_links : list) -> list:
    """
    Filter a list of PDF links by user choice.

    Params:
    - pdf_links (list): A list of dictionaries containing PDF links with keys such as 'href'.

    Returns:
    - list: A list containing the selected PDF link(s) based on the user's choice.
    """
    for index, pdf_link in enumerate(pdf_links):
        file_name = pdf_link['href'].split('/')[-1]
        print(f"{index}: {file_name}")
    
    print("Enter the index of the file name you would like to download.\n\
Enter 'a' to download all pdf files.\n\
Enter 'q' to quit.")
          
    while True:
        choice = input()

        # filter to an empty list
        if choice.lower() == 'q':
            print("Quitting...")
            return []

        # fall through without filtering
        if choice.lower() == 'a':
            break

        try:
            index = int(choice)
            if -1 < index < len(pdf_links):
                # filter to a single pdf link as per user choice
                pdf_links = [pdf_links[index]]
                print("Selected file:", (pdf_links[0])['href'].split('/')[-1])
                break
            else:
                print("Invalid choice. Please try again.")

        except ValueError:
            print("Invalid input. Please enter a number within index ranges, 'q' to quit or 'a' to download all.")

    return pdf_links 

def download_pdf(url : str, download_all : bool = False, target_pdfs_directory : str = 'pdf_extraction/pdfs') -> list:
    """
    Attempts to download PDF files from a given URL.

    Params:
    - url (str): The URL to download PDF files from.
    - download_all (bool): Whether the user will be given a choice to specify which pdfs to download, or if all will be downlaoded
    - target_pdfs_directory (str): Directory that the downloaded pdfs will be saved into

    Returns:
    - List[Tuple[str, bytes]]: A list of tuples, where each tuple contains the file name and content of a downloaded PDF file.
      If no PDF files are found or an error occurs, an empty list is returned.
    """

    print(f"Attempting to download pdf files from {url}")

    # headers simulating a browser request are necessary, without them the request status code returns 403 (forbidden)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # headers=headers prevents 403 forbidden response.status_code
    response = requests.get(url, headers=headers)

    # 200 -> good status code
    if response.status_code != 200:
        print(f"Failed to download the web page. response.status_code: {response.status_code}")
        return []

    # parses the web page for the first pdf file (assuming the molecules
    # of the month file is the first one on the web page)
    soup = BeautifulSoup(response.content, 'html.parser')
    pdf_links = soup.select('a[href$=".pdf"]')

    if not pdf_links:
        print("PDF links not found on the page.")
        return []

    # allow user to specify which pdf links to continue with
    if not download_all:
        pdf_links = filter_pdf_links_by_choice(pdf_links)

    # attempt to download pdf files and save them into a list of tuples (filename, content)
    pdf_files = []
    for pdf_link in pdf_links:
        file_name = pdf_link['href'].split('/')[-1]
        file_url = pdf_link['href']   
        pdf_response = requests.get(file_url, headers=headers)
        if pdf_response.status_code == 200:
            pdf_files.append((file_name, pdf_response.content))
            print(file_name + " downloaded successfully.")
        else:
            print("Failed to download the PDF: ", file_name)

    # save to specified directory
    if target_pdfs_directory:
        os.makedirs(target_pdfs_directory, exist_ok=True)
        for filename, pdf_content in pdf_files:
            with open(os.path.join(target_pdfs_directory, filename), 'wb') as f:
                f.write(pdf_content)
                print(f"{filename} saved into {target_pdfs_directory}.")

    return pdf_files


def main():

    # example of use
    # working url

    url = "https://drughunter.com/resource/2022-drug-approvals/"
    download_pdf(url, download_all=False, target_pdfs_directory="pdfs")

    

if __name__ == '__main__':
    main()
