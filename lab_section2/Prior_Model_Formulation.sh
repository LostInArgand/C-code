#!/bin/bash

# Compute the noncausal prediction error for the image img04.tif
echo "Computing the noncausal prediction error for the image img04.tif"
../bin/ncpe ../images/img04g.tif > output/sigmas.txt
mv ncpe_img.tif output
# mv color.tif output
echo ""


# Run this matlab script to plot the data
# matlab plot_sigma.m

