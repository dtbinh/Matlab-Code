% This script plots all the possible stations for each group of
% proxies in the Pseudoproxies folder. This will use a density plot for all
% the calibration windows

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
for window=[31,61,91]
    
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 499;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:NUM_CAL_WDW-1
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
figsopen = 1; close all;
folderlist = dir(['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/']);
folderlist([3,4,6,8,9,11,12,14,15]) = [];% Remove the ones you dont want
for i=1:(length(folderlist)-2)
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',folderlist(i+2).name];
    lat_ind = []; lon_ind=[]; num_prox = zeros(NUM_CAL_WDW,1);
    for c=1:NUM_CAL_WDW
        load([DIR_NAME, '/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/3stns_1000prox.mat']);
        [lat_ind_tmp,lon_ind_tmp]=ind2sub([90 144],indice_pool);
        lat_ind = [lat_ind; lat_ind_tmp]; lon_ind = [lon_ind; lon_ind_tmp];
        num_prox(c) = length(indice_pool);
    end
    subplot(3,1,rem(i-1,3)+1)
    values = hist3([lon_ind lat_ind],{1:length(lon) , 1:length(lat)})';
    pcolor(lon,lat,values); plotworld; colormap(flipud(hot(10))); colorbar;
    xlim([0 360]); ylim([-90 90]);
    title([strrep(DIR_NAME(58:end),'_','\_'),' - Density with all CalWdws, mean no.: ',num2str(mean(num_prox),'%.0f'),', std: ',num2str(std(num_prox),'%.0f')]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    if rem(i,3)==0 & i~=(length(folderlist)-2)
        figure;
        figsopen = figsopen + 1;
    end
                
end

for f=1:figsopen
    figure(f)
    saveas(gcf,['Plots/prox_location_density_',num2str(window),'yr_',num2str(f),'.jpg'])
    close
end

end