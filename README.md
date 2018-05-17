# Get PATRIC Functions
### Get functional annotation data from PATRIC

The purpose of this repo is to provide a set of functions that will
fetch the functional annotations and 16S sequences from PATRIC

#### Scripts

  * get_patric_data.py: Download data from PATRIC and reformat to be compatible with PICRUSt

```
usage: get_patric_data.py [-h] --output-folder OUTPUT_FOLDER [--test]

Download reference data from PATRIC

optional arguments:
  -h, --help            show this help message and exit
  --output-folder OUTPUT_FOLDER
                        Folder for downloaded genome data.
  --test                Use this flag to download a subset of the data for
                        testing.
```