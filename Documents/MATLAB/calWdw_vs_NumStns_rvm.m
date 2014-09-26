%% Correlations

figure;
for window = [31, 61, 91]

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_RVM','all_stn_corr_MRV')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
end

% Plotting EPC
subplot(3,3,1+(floor(window/30)-1)*3) % temp_corr_EPC_RV
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
ylabel(['Correlation (',num2str(window),'yrs of data)'])
if window==31 title(['RVM']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting CPS
subplot(3,3,2+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==91 xlabel('Number of Stations included in reconstruction'); end
if window==31 title(['RVM']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting MRV
subplot(3,3,3+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==31 title(['MRV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

end
suptitle([strrep(GROUP_NAME,'_','\_'),' - Ranges of Correlation percentiles'])
set(gcf, 'PaperPosition', [0 0 19 23]);
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);

%% RMSE
figure;
numstnstocompare=3:70;
for window = [31, 61, 91]

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_rmse_EPC_RV','all_stn_rmse_RVM','all_stn_rmse_MRV')
    temp_corr_EPC_RV(:,c,:) = all_stn_rmse_EPC_RV;
    temp_corr_RVM(:,c,:) = all_stn_rmse_RVM;
    temp_corr_MRV(:,c,:) = all_stn_rmse_MRV;
end

% Plotting EPC
subplot(3,3,1+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on
ylabel(['RMSE (',num2str(window),'yrs of data)'])
if window==31 title(['EPC\_RV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting CPS
subplot(3,3,2+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on
if window==91 xlabel('Number of Stations included in reconstruction'); end
if window==31 title(['RVM']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting MRV
subplot(3,3,3+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on
if window==31 title(['MRV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

end
suptitle([strrep(GROUP_NAME,'_','\_'),' - Ranges of RMSE percentiles'])
set(gcf, 'PaperPosition', [0 0 28 19]);
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);

%% Figure 4-7 _ Variance of Reconstructions
figure;
for window = [31, 61, 91]

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_var_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_var_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_var_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'var_EPC_RV','var_CPS_RV','var_MRV')
    temp_var_EPC_RV(:,c,:) = var_EPC_RV;
    temp_var_CPS_RV(:,c,:) = var_CPS_RV;
    temp_var_MRV(:,c,:) = var_MRV;
end

% Plotting EPC
subplot(3,3,1+(floor(window/30)-1)*3)
var_RV_qn = quantile(temp_var_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(var_RV_qn,1),2,size(var_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(var_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(var_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
ylabel(['Correlation (',num2str(window),'yrs of data)'])
if window==31 title(['EPC\_RV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting CPS
subplot(3,3,2+(floor(window/30)-1)*3)
var_RV_qn = quantile(temp_var_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(var_RV_qn,1),2,size(var_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(var_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(var_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==91 xlabel('Number of Stations included in reconstruction'); end
if window==31 title(['CPS\_RV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting MRV
subplot(3,3,3+(floor(window/30)-1)*3)
var_RV_qn = quantile(temp_var_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(var_RV_qn,1),2,size(var_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(var_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(var_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==31 title(['MRV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

end
suptitle([strrep(GROUP_NAME,'_','\_'),' - Ranges of Correlation percentiles'])
set(gcf, 'PaperPosition', [0 0 19 23]);
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);