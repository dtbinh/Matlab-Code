% This script will look closely at one reconstruction and its statistics

%% Setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUM_STNS = 10;
DIR_NAME = 'Pseudoproxies/glb_peof2';
RUN_NUM = 1; % Out of the thousand runs
CAL_NUM = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([DIR_NAME,'/',num2str(NUM_STNS),'stns_1000prox.mat'])
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
time = nc_varget(ts_file,'time'); % Assumes both files use the same time
ts = nc_varget(ts_file,'ts')-273.15; % To Celsius
pr = nc_varget(pr_file,'pr');
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nE] = min(abs(lon-240));
[~,nW] = min(abs(lon-190));

% Annual (Jul-Jun) Means and Anomalies
jul_jun_fmt = 7:5994;
ts=ts(jul_jun_fmt,:,:);
pr=pr(jul_jun_fmt,:,:);
time=time(jul_jun_fmt);

ats=zeros(size(ts,1)/12,size(ts,2),size(ts,3));
apr=zeros(size(pr,1)/12,size(pr,2),size(pr,3));
for y=1:length(time)/12
    ats(y,:,:)=mean(ts((12*(y-1)+1):(y*12),:,:),1);
    apr(y,:,:)=mean(pr((12*(y-1)+1):(y*12),:,:),1);
end

trend = zeros(2,size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        trend(:,i,j) = regress(ats(:,i,j), [ones(499,1) [1:length(time)/12]']);
        ats(:,i,j) = ats(:,i,j) - [1:length(time)/12]'*trend(2,i,j);
    end
end
for i=1:size(apr,2)
    for j=1:size(apr,3)
        trend(:,i,j) = regress(apr(:,i,j), [ones(499,1) [1:length(time)/12]']);
        apr(:,i,j) = apr(:,i,j) - [1:length(time)/12]'*trend(2,i,j);
    end
end

ats_anmn = mean(ats,1);
apr_anmn = mean(apr,1);
for y=1:length(time)/12
    ats(y,:,:)=ats(y,:,:)-ats_anmn;
    apr(y,:,:)=apr(y,:,:)-apr_anmn;
end
ats_anmn=squeeze(ats_anmn);
apr_anmn=squeeze(apr_anmn);
n34_ind = mean(mean(ats(:,nS:nN,nW:nE),3),2);

clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

VAR_WDW = 30; % Moving window for moving variance is 30 Years
window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];
CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];
load([DIR_NAME,'/tonsofstats.mat'],'CAL_WDW','all_stn_CPS_RV','all_stn_EPC_RV','all_stn_MRV',...
                                   'all_stn_CPS','all_stn_EPC',...
                                   'all_stn_corr_CPS','all_stn_rmse_CPS','all_stn_corr_EPC','all_stn_rmse_EPC',...
                                   'all_stn_corr_CPS_RV', 'all_stn_corr_EPC_RV', 'all_stn_corr_MRV',...
                                   'all_stn_rmse_CPS_RV', 'all_stn_rmse_EPC_RV', 'all_stn_rmse_MRV');
%% Plotting Proxy locations
clf;
subplot(3,1,1)
plotworld;
hold on; scatter(lon(stn_lon(RUN_NUM,:)),lat(stn_lat(RUN_NUM,:)),'go','filled'); hold off;
xlim([0 360]); ylim([-90 90]);
title(['Proxy locations for run ',num2str(RUN_NUM),' of ',strrep(DIR_NAME(15:end),'_','\_'),', NUM\_STNS = ',num2str(NUM_STNS)])

% Plotting Time series of recontructions

subplot(3,1,2); hold on;
plot(n34_ind_RV,'k','LineWidth',2);
plot(squeeze(all_stn_EPC_RV(NUM_STNS,CAL_NUM,RUN_NUM,:)),'r');
plot(squeeze(all_stn_CPS_RV(NUM_STNS,CAL_NUM,RUN_NUM,:)),'b');
plot(squeeze(all_stn_MRV(NUM_STNS,RUN_NUM,:)),'g'); 
hold off;
legend(['N34\_ind\_RV'], ...
       ['stn\_EPC\_RV, Corr: ',num2str(all_stn_corr_EPC_RV(NUM_STNS,CAL_NUM,RUN_NUM),3),', RMSE: ',num2str(all_stn_rmse_EPC_RV(NUM_STNS,CAL_NUM,RUN_NUM),3)],...
       ['stn\_CPS\_RV, Corr: ',num2str(all_stn_corr_CPS_RV(NUM_STNS,CAL_NUM,RUN_NUM),3),', RMSE: ',num2str(all_stn_rmse_CPS_RV(NUM_STNS,CAL_NUM,RUN_NUM),3)],...
       ['stn\_MRV, Corr: ',num2str(all_stn_corr_MRV(NUM_STNS,RUN_NUM),3),', RMSE: ',num2str(all_stn_rmse_MRV(NUM_STNS,RUN_NUM),3)],...
       'location','northeast');
title(['Running Variance timeseries of run: ',num2str(RUN_NUM),' of ',strrep(DIR_NAME(15:end),'_','\_'), ...
       ', NUM\_STNS = ',num2str(NUM_STNS),', CAL\_WDW = ',num2str(CAL_WDW(CAL_NUM,1)),':',num2str(CAL_WDW(CAL_NUM,end))])

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
ylim([0.3 2.5]); grid on;
% Plotting time series
subplot(3,1,3);
hold on;
plot(n34_ind,'k','LineWidth',2.5);
plot(squeeze(all_stn_EPC(NUM_STNS,CAL_NUM,RUN_NUM,:)),'r','LineWidth',2);
plot(squeeze(all_stn_CPS(NUM_STNS,CAL_NUM,RUN_NUM,:)),'b','LineWidth',1);
hold off;
title(['Timeseries of run: ',num2str(RUN_NUM),' of ',strrep(DIR_NAME(15:end),'_','\_'), ...
       ', NUM\_STNS = ',num2str(NUM_STNS),', CAL\_WDW = ',num2str(CAL_WDW(CAL_NUM,1)),':',num2str(CAL_WDW(CAL_NUM,end))])
legend(['N34\_ind'], ...
       ['stn\_EPC, Corr: ',num2str(all_stn_corr_EPC(NUM_STNS,CAL_NUM,RUN_NUM),3),', RMSE: ',num2str(all_stn_rmse_EPC(NUM_STNS,CAL_NUM,RUN_NUM),3)],...
       ['stn\_CPS, Corr: ',num2str(all_stn_corr_CPS(NUM_STNS,CAL_NUM,RUN_NUM),3),', RMSE: ',num2str(all_stn_rmse_CPS(NUM_STNS,CAL_NUM,RUN_NUM),3)],...
       'location','northeast');
ylim([-3 4]); grid on;
saveas(gcf,['Plots/examSpecifRecons_',DIR_NAME(15:end),'_nStn',num2str(NUM_STNS),'_rn',num2str(RUN_NUM),'_clWdw',num2str(CAL_WDW(CAL_NUM,1)),':',num2str(CAL_WDW(CAL_NUM,end)),'.jpg']);