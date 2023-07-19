import requests
from bs4 import BeautifulSoup

def download_pdf(url:str):
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

def main():
    # example of use
    # working url

    url = "https://drughunter.com/molecules-of-the-month/2023/may-2023/"
    pdf_list = download_pdf(url)
    
    print(f"Extract {len(pdf_list)} files in total.")
    for filename, pdf in pdf_list:
        print(f"{filename} extracted.")
    

if __name__ == '__main__':
    main()
