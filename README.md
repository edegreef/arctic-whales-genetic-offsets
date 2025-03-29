![github_header_arcticwhales_offset](https://github.com/user-attachments/assets/90509c18-e4dd-44fc-93ab-bc196bf24f54)

This is a repository for scripts used in analyzing Canadian Arctic bowhead whale (*Balaena mysticetus*), narwhal (*Monodon monoceros*), and beluga whale (*Delphinapterus leucas*) resequencing data for analyses to run gradient forest models and calculate genetic offsets. These scripts were executed on computing resources through the University of Manitoba and the Digitial Research Alliance of Canada. 

A complementary repository "Arctic whales resequencing" provides scripts for starting with raw sequence files (corresponding to https://doi.org/10.1111/gcb.17528) and can be found here: https://github.com/edegreef/arctic-whales-resequencing. 
SNP data and corresponding metadata used in the gradient forest and genetic offset analyses are available at https://doi.org/10.5281/zenodo.15105821.

Scripts in this current repository include:
* Creating map with combined species range shapefiles
* SNP filtering pipeline
* Imputing with *beagle* and thinning dataset with *vcftools*
* Extracting climate data with *sdmpredictors*
* Running genoome scans and lfmms in *LEA*
* SNP prep and formating for Gradient Forest
* Pipeline for gradient forest analyses with *gradientForest*
* Scaling and combining genetic offsets
