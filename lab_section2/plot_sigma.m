
clear all
close all
data = load("/Users/pradithaalwis/Documents/Academics/Purdue - PhD/Semester 1/ECE 64100/Labs/MAP_Image_Restoration/C-code/lab_section2/output/sigmas.txt");
data = data';

p = data(1,:);
sigma = data(2,:);

plot( p, sigma );
title("Plot of ML estimate of the scalar parameter \sigma vs. p");
xlabel('p');
ylabel('\sigma');

