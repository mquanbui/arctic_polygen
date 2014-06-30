clear all
close all
clc


centers = importdata('centers.txt');

radius = 4;

hold on
for i = 2:size(centers,1)
    circle(centers(i,2:3),radius,1000,'-');
    scatter(centers(i,2),centers(i,3), 'black', 'fill')
end
