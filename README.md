# DrugHunter Molecule Extractor

The DrugHunter website publishes many high quality datasets. Specifically
the yearly [Drug Approvals](https://drughunter.com/resource_category/approved-drug-reviews/) and monthly [Molecules of the Month](https://drughunter.com/molecules-of-the-month/). These sets are well-curated and provide a useful source of validated molecules. Unfortunately
all the molecules are presented strictly contained within pdf pages - making their extraction into computer-readable format a non-trivial taks.

This repo provides the tools that allow the extraction of these molecules from the webpage - either from a provided url,
or through a workflow that extracts all Molecules of the Month molecules within a given year. All extracted information is exported into
a timestamped 

The script can also easily be used to extract molecules within pdfs from any provided url. See Usage.

## Installation

Clone repository and install all the dependencies
```bash
git clone https://github.com/deimos1078/drughunter-molecule-extractor
cd drughunter-molecule-extractor
pip install -r requirments.txt
```
Install modified Molscribe
```bash
git clone https://github.com/deimos1078/MolScribe
cd MolScribe
python setup.py install
cd ..
```

Install modified decimer_segmentation
```bash
git clone https://github.com/deimos1078/DECIMER-Image-Segmentation
cd DECIMER-Image-Segmentation
pip install .
cd ..
```

Feel free to remove the DECIMER and MolScribe repos once they're installed.

If all goes well, Drug Hunter extractor should now be usable

### Note about rdkit version
It's highly recommended to use the version of rdkit specified within requirments
Both decimer and molscribe use rdkit to generate smiles and using older versions
may lead to unexpected effects in the resulting smiles generation (unwanted hybridization for example)
## Usage

Invoking help:
```bash
python3 drughunter_extractor.py -h
```
will produce:
```text
usage: drughunter_extractor.py [-h] [-y YEAR] [-m MONTH] [-u URL] [--seg_dir SEG_DIR] [--decimer_off] [--text] [--direction DIRECTION] [--separator SEPARATOR]

DrugHunter extractor

options:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  (int) targeted year of drughunter molecules of the month set
  -m MONTH, --month MONTH
                        (str) targeted month range of the molecules of the month set, input either two numbers separated by a dash or a single number (borders of the range are included)
  -u URL, --url URL     (str) url of webpage with targeted set (in case the format of drughunter url changes, which is likely)
  --seg_dir SEG_DIR     (str) directory that the segmented segments will be saved into, if unspecified, segments will not be saved
  --decimer_off         Turns off decimer complementation
  --text                Turns on text extraction
  --direction DIRECTION
                        Specifies in which direction the text is from the molecules
  --separator SEPARATOR
                        Specifies which separator is used in the document to separate name and target.
```

### Extract from url (without text)

-u, --url  Specify the url containing pdf (or pdfs) that the script will attempt to extract molecules from 

```bash
python3 drughunter_extractor.py --url https://drughunter.com/resource/2022-drug-approvals/
```

The script will attempt to access the webpage and then proceed to list all links to pdf files on the site
You can either select a specific pdf link or download all of them for the proceeding extraction


```text
Attempting to download pdf files from https://drughunter.com/resource/2022-drug-approvals/
0: DH-2022-Small-Drug-Approvals-v3.pdf
1: DH-2022-Large-Drug-Approvals-_R1.pdf
2: EC-edits_DH-2022-First-in-Class-Small-Molecules-R1..pdf
3: DH-2022-First-in-Class-Large-Molecules.pdf
Enter the index of the file name you would like to download.
Enter 'a' to download all pdf files.
Enter 'q' to quit.
```
 
Input:
```bash
2
```

```text
Selected file: EC-edits_DH-2022-First-in-Class-Small-Molecules-R1..pdf
EC-edits_DH-2022-First-in-Class-Small-Molecules-R1..pdf downloaded successfully.
```

Now simply wait for the script to perform segmentation, recognition and validation.
The results will be exported into a csv in the results directory.


### Extract from url (with text) (recommended only for extracting text out of DrugHunter)

use --direction with either "down" or "right" to specify where the text you're looking to extract is relative to the molecule

use --separator if the text information contains a "NAME SEPARATOR TARGET" line to specify the separator

```bash
python3 drughunter_extractor.py --url https://drughunter.com/resource/2022-drug-approvals/ --text --direction "down" --separator ','
```


### Extract all Molecules of the Month within a given year

-y, --year Specify the year that the Molecules of the Month sets you're targeting were published in

Use the --text switch to specify whether you'd like the script to attempt to extract text information form the pdfs as well

```bash
python3 drughunter_extractor.py --year 2023 --text
```

Please note that at the time of writing both Decimer and MolScribe generate quite a bit of warnings.
This is expected behaviour so as long as the script exports results, everything is working as intended.

### Extract Molcules of the Month within a specified month range

-m, --month Specify the range of months that the Molecules of the Month sets you're targeting were published in


Ex. I want to download all sets published between February(2) and September(9) of 2022, but no text
```bash
python3 drughunter_extractor.py --year 2022 --month 2-9
```

Ex. I want to download the May 2023 set, with text
```bash
python3 drughunter_extractor.py --year 2022 --month 5 --text
```

Ex. I want to download the June 2023 set, without using decimer to complement the results
```bash
python3 drughunter_extractor.py --year 2022 --month 6 --decimer_off
```

## Workflow and used libraries

1) The webpage is accessed through [requests](https://pypi.org/project/requests/)
2) [BeautifulSoup](https://pypi.org/project/beautifulsoup4/) is used to gather all pdf links on the webpage
3) Pdfs are downloaded using those links
4a) Pdfs are segmented into individual images using [decimer-image-segmentation](https://github.com/Kohulan/DECIMER-Image-Segmentation/tree/master)
4b) If the molecules of the months set is targeted, segmentation is done using rectangle detection instead
5) If --text option is on, the bounding boxes of these segments are used to attempt text extraction as well (using fitz from [pymupdf](https://pymupdf.readthedocs.io/en/latest/module.html))
6) The segments are recognized using [MolScribe](https://github.com/thomas0809/MolScribe)
7) Inchikeys are gathered from the recognized smiles using [rdkit](https://www.rdkit.org/)
8) Inchikeys are searched for in the [Unichem](https://www.ebi.ac.uk/unichem/) database for connectivity in order to access their validity using [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client)
9) Segments that were not validated by Unichem are recognized by [Decimer-Image_Transformer](https://github.com/Kohulan/DECIMER-Image_Transformer). This is done because Decimer and MolScribe are good at recognizing different molecules, and because decimer is significantly slower, it is prefferable
to use it on only a neccessary portion of the segments
10) Validation again
11) The results are filtered so that no duplicate inchikeys are present
12) The results are exported into a csv file using [pandas dataframe](https://pandas.pydata.org/)

