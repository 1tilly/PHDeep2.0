## Current state
This code has not been tested since the refactoring. This is going to happen, when a large chunk of the code has been moved. Further, there are a few placeholder files/direcotires present, as the first step in the refactoring was the consideration of modules and components in this repository.

# PHDeep - Pulmonary Hypertension Deep Learning pipeline
This project is a refactored version of my PhD project, based on the thesis "Deep Learning of regulatory sequence variation in Pulmonary Arterial Hypertension". It is still strongly under development.

# Project structure
This project is modularised to be as adjustable as possible, while still (hopefully) easy to understand. Most module names should be self-explanatory. Something to point out are prediction and post_prediction. The former holding functions to actually run across either the reference genome or through vcf files, fetching the reference around the variants on the go. For this, the data_loading module is used to load the files. 
Post_prediction is for the handling of predicted epigenetic marks. Here, you will find the functions that compute the scores and make them accessible to the statistical/association testing. 

# Future work
These modules will be connected in Nextflow Pipelines, so that it is possible to run the pipeline from training-set generation over model training to prediction and post_prediction. 
