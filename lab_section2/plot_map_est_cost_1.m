
clear all
close all
path = "/Users/pradithaalwis/Documents/Academics/Purdue - PhD/Semester 1/ECE 64100/Labs/MAP_Image_Restoration/C-code/lab_section2/output/cost_MAP_est_blurred_noisy.txt";
data = load(path);
data = data';

k = data(1,:);
cost = data(2,:);

plot( k, cost );
title("Plot of cost of MAP estimate of X vs the iteration number");
xlabel('Iteration Number');
ylabel('Cost');

