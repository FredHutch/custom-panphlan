# custom-panphlan
### Make a custom PanPhlAn database

The purpose of this repo is to provide a set of functions that will
make it easier to build a PanPhlAn database that is specific to one
particular experiment.

#### Workflow

##### Database Build Time

  * Download data from PATRIC
  * Reformat to be compatible with PanPhlAn

##### Experiment Analysis Time

  * Recruit 16S sequences from PATRIC 16S references (and others)
  * Build a tree, including all sequence variants
  * Make a custom PanPhlAn database
  * Predict functional capacity for each sample in the experiment
