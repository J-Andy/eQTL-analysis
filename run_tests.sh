#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p $(pwd)/output

# Run the Singularity container and execute the unit tests
singularity exec --bind $(pwd)/output:/usr/src/app/output r-analysis.sif Rscript /usr/src/app/unit_tests.R