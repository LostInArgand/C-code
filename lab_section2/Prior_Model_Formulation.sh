#!/bin/bash

# Compute the noncausal prediction error for the image img04.tif
echo "Computing the noncausal prediction error for the image img04.tif"
../bin/ncpe ../images/img04g.tif > output/sigmas.txt
mv ncpe_img.tif output
# mv color.tif output
echo ""

# Perform ICD optimization on the noisy image img04.tif
echo "Performing ICD optimization on the noisy image img04.tif"
../bin/ICD_opt ../images/img04g.tif > output/logs.txt
mv noisy_img.tif output
mv MAP_est_img_sigma_x_1.tif output
mv MAP_est_img_sigma_x_5.tif output
mv MAP_est_img_sigma_x_1_5.tif output

mv costs_sigma_1.txt output
mv costs_sigma_5.txt output
mv costs_sigma_1_5.txt output

# Perform ICD optimization on the blurred and noisy image of img04.tif
echo "Perform ICD optimization on the blurred and noisy image of img04.tif"
../bin/ICD_BN_opt ../images/img04g.tif > output/cost_MAP_est_blurred_noisy.txt
mv noisy_blurred_img.tif output
mv MAP_est_blurred_noisy_img.tif output

# Perform ICD optimization with a non guassian prior on the blurred and noisy image of img04.tif
echo "Perform ICD optimization with a non gaussian prior on the blurred and noisy image of img04.tif"
../bin/non_gauss ../lab_section2/output/noisy_blurred_img.tif > output/cost_MAP_est_non_gauss_blurred_noisy.txt
mv MAP_est_non_gauss_img.tif output

# Perform ICD optimization with a non guassian prior with majorization on the blurred and noisy image of img04.tif
echo "Perform ICD optimization with a non gaussian prior with majorization on the blurred and noisy image of img04.tif"
../bin/non_gauss_with_maj ../lab_section2/output/noisy_blurred_img.tif > output/cost_MAP_est_non_gauss_with_maj_blurred_noisy.txt
mv MAP_est_non_gauss_with_maj_img.tif output


# Run this matlab script to plot the data
/Applications/MATLAB_R2024a.app/bin ./output/plot_sigma.m

