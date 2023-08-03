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

## Usage
The drughunter_extractor.py has two modes

### Extract from url

-u, --url  Specify the url containing pdf (or pdfs) that the script will attempt to extract molecules from 

```bash
python3 drughunter_extractor.py --url https://drughunter.com/resource/2022-drug-approvals/
```

The script will attempt to access the webpage and then proceed to list all links to pdf files on the site
You can either select a specific pdf link or download all of them for the proceeding extraction


```bash
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
$ 2
```

```bash
Selected file: EC-edits_DH-2022-First-in-Class-Small-Molecules-R1..pdf
EC-edits_DH-2022-First-in-Class-Small-Molecules-R1..pdf downloaded successfully.
```

Now simply wait for the script to perform segmentation, recognition and validation.
The results will be exported into a csv in the results directory.

### Extract all Molecules of the Month within a given year

-y, --year Specify the year that the Molecules of the Month sets you're targeting were published in

Just specify a year and the script will do the rest (coffee time)

```bash
python3 drughunter_extractor.py --year 2023
```

Please note that at the time of writing both Decimer and MolScribe generate quite a bit of warnings.
This is expected behaviour so as long as the script exports results, everything is working as intended.

## Workflow and used libraries

1) The webpage is accessed through [requests](https://pypi.org/project/requests/)
2) [BeautifulSoup](https://pypi.org/project/beautifulsoup4/) is used to gather all pdf links on the webpage
3) Pdfs are downloaded using those links
4) Pdfs are segmented into individual images using [decimer-image-segmentation](https://github.com/Kohulan/DECIMER-Image-Segmentation/tree/master)
5) The segments are recognized using [MolScribe](https://github.com/thomas0809/MolScribe)
6) Inchikeys are gathered form the recognized smiles using [rdkit](https://www.rdkit.org/)
7) Inchikeys are searched for in the [Unichem](https://www.ebi.ac.uk/unichem/) database for connectivity in order to access their validity using [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client)
8) Segments that were not validated by Unichem are recognized by [Decimer-Image_Transformer](https://github.com/Kohulan/DECIMER-Image_Transformer). This is done because Decimer and MolScribe are good at recognizing different molecules, and because decimer is significantly slower, it is prefferable
to use it on only a neccessary portion of the segments
9) Validation again
10) The results are exported into a csv file using [pandas dataframe](https://pandas.pydata.org/), containing the source pdf file name, recognized smiles, inchi, inchikey and whether the inchikey was validated or not

