# DrugHunter Molecule Extractor

The DrugHunter website publishes many high quality datasets. Specifically
the yearly Drug Approvals (https://drughunter.com/resource_category/approved-drug-reviews/) and monthly Molecules of the Month (https://drughunter.com/molecules-of-the-month/). These sets are well-curated and provide a useful source of validated molecules. Unfortunately
all the molecules are presented strictly contained within pdf pages - making their extraction into computer-readable format a non-trivial taks.

This library provides the tools that allow the extraction of these molecules from the webpage - either from a provided url,
or through a workflow that extracts all Molecules of the Month molecules within a given year.


