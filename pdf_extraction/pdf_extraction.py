import requests
from bs4 import BeautifulSoup

def filter_pdf_links_by_choice(pdf_links:list):
    # print out all filenames to choose from
    for index, pdf_link in enumerate(pdf_links):
        file_name = pdf_link['href'].split('/')[-1]
        print(f"{index}: {file_name}")
    
    print("Enter the index of the file name you would like to download.\n\
Enter 'a' to download all pdf files.\n\
Enter 'q' to quit.")
          
    while True:
        choice = input()

        if choice.lower() == 'q':
            print("Quitting...")
            return []

        if choice.lower() == 'a':
            break

        try:
            index = int(choice)
            if -1 < index < len(pdf_links):
                pdf_links = [pdf_links[index]]
                print("Selected file:", (pdf_links[0])['href'].split('/')[-1])
                break
            else:
                print("Invalid choice. Please try again.")

        except ValueError:
            print("Invalid input. Please enter a number within index ranges, 'q' to quit or 'a' to download all.")

    return pdf_links 

def download_pdf(url:str, download_all:bool = False):
    """
    Attempt to download PDF files from a given URL.

    Params:
    - url (str): The URL to download PDF files from.

    Returns:
    - List[Tuple[str, bytes]]: A list of tuples, where each tuple contains the file name and content of a downloaded PDF file.
      If no PDF files are found or an error occurs, an empty list is returned.

    The function sends a GET request to the specified URL with headers simulating a browser request to avoid a 403 (forbidden) status code.
    If the response status code is 200, it parses the web page using BeautifulSoup to find all PDF file links.
    For each PDF link found, it retrieves the file name and URL, and sends another GET request to download the file.
    If the download is successful (status code 200), it appends a tuple containing the file name and file content to the `pdf_files` list.
    If any errors occur during the download process, error messages are printed.

    Finally, the function returns the list of downloaded PDF files.
    """

    print(f"Attempting to download pdf files from {url}")

    # headers simulating a browser request are necessary, without them the request status code returns 403 (forbidden)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # headers=headers prevents 403 forbidden response.status_code
    response = requests.get(url, headers=headers)

    if response.status_code != 200:
        print("Failed to download the web page.")
        return []

    # parses the web page for the first pdf file (assuming the molecules
    # of the month file is the first one on the web page)
    soup = BeautifulSoup(response.content, 'html.parser')
    pdf_links = soup.select('a[href$=".pdf"]')

    if not pdf_links:
        print("PDF links not found on the page.")
        return []

    if not download_all:
        pdf_links = filter_pdf_links_by_choice(pdf_links)

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
    return pdf_files


def main():
    from pathlib import Path
    import os
    # example of use
    # working url

    url = "https://drughunter.com/resource/2022-drug-approvals/"
    pdf_list = download_pdf(url, download_all=True)
    
    os.makedirs("pdf_extraction/pdfs", exist_ok=True)

    print(f"Extract {len(pdf_list)} files in total.")
    for filename, pdf_content in pdf_list:
        with open(os.path.join("pdf_extraction/pdfs", filename), 'wb') as f:
            f.write(pdf_content)
        print(f"{filename} extracted.")

    

if __name__ == '__main__':
    main()
