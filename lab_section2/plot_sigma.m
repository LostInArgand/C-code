
clear all
close all
load output/sigmas.txt
data = data';

p = data(1,:);
sigma = data(2,:);

plot( p, sigma );
xlabel('p');
ylabel('sigma');

