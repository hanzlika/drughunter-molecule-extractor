# drughunter molecule extractor
 A miniproject focused on creating a script for repeated extraction of molecules and attached text information from pdf files published on https://drughunter.com/, primarily on the sets published on  https://drughunter.com/molecules-of-the-month/ and  https://drughunter.com/resource_category/approved-drug-reviews/ using OCSR.

 The workflow should consist of these steps:
	1) PDF extraction
	2) PDF segmentation
	3) Optical Chemical Structure Recognition (OCSR, most likely with DECIMER)
	4) Validation through cross-reference with PubChem or UniChem
	5) Results export as .csv