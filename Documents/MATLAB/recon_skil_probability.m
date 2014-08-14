% This script will plot the 'tonsofstats.mat' file data as ranges of the
% calibration windows. The file data were the ones produced by Katana, and
% are stored on the server

% This script needs mexcdf to be installed in matlab to run, the
% function 'movingvar', 'plotworld', and the folder DataFiles

ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
numstnstocompare = 70; NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
% Calibration windows set to being 10 overlapping windows over 499 years
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
%% Plotting
GROUP_NAME = 'glb_ts';
for window = [31,61,91]
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])
    
    for group_size=3:numstnstocompare
        
        all_nstat_yrs = nan([10,NUM_TRIALS]);
        all_nstat_nstns = nan([10,NUM_TRIALS]);
        corrs = nan([10,NUM_TRIALS]);
        for c=1:size(CAL_WDW,1)
            DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'stn_lat','stn_lon','indice_pool');

            % Compare the lat and lon with nonstat map, save no. of nonstat, save
            % nonstat year nums, and no. of nonstat points in reconstruction Plot
            % numbers of non-stat stations/av years against Reconstruction skill. Also
            % find percentages of the runs that have the numbers of the non-stat
            % prox/av yrs

            nstat_yrs = nan(size(stn_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs(tr,m) = nonstat_tsmap(stn_lat(tr,m),stn_lon(tr,m));
                end
            end

            nstat_avyrs = mean(nstat_yrs,2);
            nstat_nstns = sum(nstat_yrs > ceil(0.1*(NUM_YRS-window)),2)/group_size;
            all_nstat_yrs(c,:) = nstat_avyrs;
            all_nstat_nstns(c,:) = nstat_nstns;
            
            
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
 'all_stn_corr_CPS_RV', 'all_stn_corr_EPC_RV', 'all_stn_corr_CPS','all_stn_corr_EPC','all_stn_corr_MRV');
            corrs(c,:) = squeeze(all_stn_corr_CPS_RV(group_size,:))';
        end
% Plots
% NUM_BINS = 10; 
% subplot(2,1,1); hist(nstat_avyrs,NUM_BINS); xlim([0 200]); title('Av number of nstat years per reconstruction');
% subplot(2,1,2); hist(nstat_nstns,NUM_BINS); xlim([0 1]); title('Percentage of nstat stns per reconstruction');

data = [all_nstat_nstns(:), corrs(:)]; 
hist3(data); xlabel('Percentage of nstat stns'); ylabel('Correlation of reconstructions');
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        
    end
end

%% Other Plotting

subplot(2,1,2);
hist(all_nstat_nstns',[0:0.025:1.0]+0.025/2); grid minor
xlim([0 0.5]); xlabel(' Percentage of non-stationary stations ');
% h = findobj(gca,'Type','Patch');
% set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
ylabel('Number of reconstructions (out of 1000)')
title('Percentage of nstat stns per reconstruction');
subplot(2,1,1);
hist(all_nstat_yrs',[0:5:500]+5/2); xlim([0 100]); grid minor;
% h = findobj(gca,'Type','Patch');
% set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
xlabel('Average number of non-stationary years per proxy');
ylabel('Number of reconstructions (out of 1000)')
title('Av number of nstat years per proxy used in the reconstructions');

set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 28 19]); %x_width=19cm y_width=28cm
suptitle([strrep(DIR_NAME(58:end),'_','\_'),'  Non-stationary pseudoproxies in reconstructions, group\_size=',num2str(group_size)]);
saveas(gcf,['Plots/hist(nstat_nstns_nstat_avyrs)',num2str(window),'rcorwdw_grp',num2str(group_size),'_ts.jpg'])