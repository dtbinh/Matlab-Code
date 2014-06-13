% This script plots the world map over the figure in focus
% It requires the lat and lon to be loaded in the plot, otherwise weird
% results will appear

load('DataFiles/coast_v2.mat');
hold on;
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
hold off
shading flat
xlabel('Longitude {\circ}E','FontSize',14);
ylabel('Latitude {\circ}N','FontSize',14);
clear lati_e lati_w long_e long_w