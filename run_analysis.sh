#!/bin/bash

# Create an output directory if it doesn't exist
mkdir -p $(pwd)/output

# Run the container and bind the output directory
singularity run --bind $(pwd)/output:/usr/src/app/output r-analysis.sif
