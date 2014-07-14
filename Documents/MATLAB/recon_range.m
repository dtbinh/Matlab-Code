% This script will plot the 'tonsofstats.mat' file data as ranges of the
% calibration windows. The file data were the ones produced by Katana, and
% are stored on the server

% This script needs mexcdf to be installed in matlab to run, the
% function 'movingvar', 'plotworld', and the folder DataFiles

ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
numstnstocompare = 70; NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;

%% Plotting
window = 31;
GROUP_NAME = 'glb_ts_nstat_sigpcd';
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

%% Series 1 Proxies plotting

%% Series 2 Proxies Plotting

CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
end

% Plotting EPC
subplot(1,3,1)
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
ylabel('Correlation')
title(['EPC\_RV'])
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting CPS
subplot(1,3,2)
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
xlabel('Number of Stations included in reconstruction');
title(['CPS\_RV'])
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting MRV
subplot(1,3,3)
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast')
xlim([0,70]); ylim([0,1]); grid on
title(['MRV'])
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

suptitle([strrep(GROUP_NAME,'_','\_'),' & ',num2str(window),'yrwdw - Ranges of Correlation percentiles'])
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast');
set(gcf, 'PaperPosition', [0 0 28 19]);
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

        
set(legendH, 'FontSize',10);
saveas(gcf,['Plots/calWdwCorr_vs_NumStns_',GROUP_NAME,num2str(window),'yr.jpg'])
