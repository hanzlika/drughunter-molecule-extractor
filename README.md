# drughunter molecule extractor
 A miniproject focused on creating a script for repeated extraction of molecules and attached text information from pdf files published on https://drughunter.com/, primarily on the sets published on  https://drughunter.com/molecules-of-the-month/ and  https://drughunter.com/resource_category/approved-drug-reviews/ using OCSR.

 The workflow should consist of these steps: <br />
	1) PDF extraction <br />
	2) PDF segmentation using https://github.com/Kohulan/DECIMER-Image-Segmentation/tree/master <br />
	3) Optical Chemical Structure Recognition (OCSR, most likely with DECIMER) <br />
	4) Validation through cross-reference with PubChem or UniChem <br />
	5) Results export as .csv <br />