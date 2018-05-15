# custom-piecrust
### Make a custom PICRUSt database

The purpose of this repo is to provide a set of functions that will
make it easier to build a PICRUSt database that is specific to one
particular experiment.

#### Workflow

##### Database Build Time

  * Download data from PATRIC
  * Reformat to be compatible with PICRUSt

##### Experiment Analysis Time

  * Recruit 16S sequences from PATRIC 16S references (and others)
  * Build a tree, including all sequence variants
  * Make a custom PICRUSt database
  * Predict functional capacity for each sample in the experiment


#### Scripts

  * get_patric_data.py: Download data from PATRIC and reformat to be compatible with PICRUSt
  * make_picrust_db.py: Given a Newick tree and the pathway annotation from PATRIC, make a custom PICRUSt database
