% This script plots all the possible stations for each group of
% proxies in the Pseudoproxies folder.

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
window = 31;
figsopen = 1; close all;
folderlist = dir(['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/']);
for i=1:(length(folderlist)-2)
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',folderlist(i+2).name];
    
    load([DIR_NAME, '/CalWdw:1-',num2str(window),'/50stns_1000prox.mat']);
    subplot(3,1,rem(i-1,3)+1)
    plotworld; hold on; scatter(lon(stn_lon(:)),lat(stn_lat(:)),'g.'); hold off;
    xlim([0 360]); ylim([-90 90]);
    title([strrep(DIR_NAME(58:end),'_','\_'),' CalWdw:1-',num2str(window)]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    if rem(i,3)==0 & i~=(length(folderlist)-2)
        figure;
        figsopen = figsopen + 1;
    end
                
end

for f=1:figsopen
    figure(f)
    saveas(gcf,['Plots/prox_location_',num2str(window),'yr_',num2str(f),'.jpg'])
    close
end
