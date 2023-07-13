import requests
from bs4 import BeautifulSoup
import urllib.request
import os

def download_pdf(target_year, target_month, target_directory):
    """Attempt to download a pdf file from drughunter molecules
    of the month based on given target year and month
    
    As of 13/7/2023 as seen in 13-7-2023-test.log
    succesfully downloaded all available molecules of the month pdfs
    (published from february 2020 to may 2023)
    """

    # general target format
    url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{target_month}-{target_year}"

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
        pdf_link = soup.select_one('a[href$=".pdf"]')

        if pdf_link:
            file_name = pdf_link['href'].split('/')[-1]
            file_url = pdf_link['href']
            file_path = os.path.join(target_directory, file_name)

            with open(file_path, 'wb') as file:
                # once again necessary to include the headers=headers file, otherwise drughunter will forbid access
                pdf_response = requests.get(file_url, headers=headers)
                file.write(pdf_response.content)

            print(file_name + " downloaded successfully.")
        else:
            print("PDF link not found on the page.")
    else:
        print("Failed to download the PDF.")

def main():
    target_dir = 'pdf_extraction\pdfs'
    for target_year in range(2020, 2024):
        for target_month in ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december']:
            print(f"Attempting to download {target_month}-{target_year}")

            download_pdf(target_year, target_month, target_dir)

if __name__ == '__main__':
    main()
