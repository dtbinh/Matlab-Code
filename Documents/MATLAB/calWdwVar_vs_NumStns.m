% This script will calculate the range of the error in the calibration windows 
% for varying station numbers, using staggered windows. It also outputs a
% mat file containing alot of statistics produced for each method,
% num_stns, and calibration window.

% Should look like 10 percentile lines for each calibration window, vs the num_stns.

% This script needs mexcdf to be installed in matlab to run, the
% function 'movingvar', 'plotworld', and the folder DataFiles

%% Setup
clear; clc;
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
window = 61; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];

%% Loading proxies
GROUP_NAME = 'glb_ts';
% DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
numstnstocompare = [3:70];                                     %%%%%%%%
NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
% Calibration windows set to being 10 overlapping windows over 499 years

%% Beginning the Loop
for window = [31, 61, 91]
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
    NUM_CAL_WDW = 10; clear CAL_WDW;
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
    end
    for c=1:size(CAL_WDW,1)

    var_MRV=nan(max(numstnstocompare), NUM_TRIALS,'single');
    var_CPS_RV=nan(max(numstnstocompare), NUM_TRIALS,'single');
    var_EPC_RV=nan(max(numstnstocompare), NUM_TRIALS,'single');
    var_RVM=nan(max(numstnstocompare), NUM_TRIALS,'single');
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'])

    for g = numstnstocompare
        for i=1:NUM_TRIALS
            var_MRV(g,i) = var(squeeze(all_stn_MRV(g,i,:)));
            var_CPS_RV(g,i) = var(squeeze(all_stn_CPS_RV(g,i,:)));
            var_EPC_RV(g,i) = var(squeeze(all_stn_EPC_RV(g,i,:)));
            var_RVM(g,i) = var(squeeze(all_stn_RVM(g,i,:)));
        end
    end
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
    save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
         'all_stn_MRV','all_stn_corr_MRV','all_stn_rmse_MRV', ...
         'all_stn_RVM','all_stn_corr_RVM','all_stn_rmse_RVM', ...
         'DIR_NAME','CAL_WDW','all_stn_EPC','all_stn_corr_EPC','all_stn_rmse_EPC','all_stn_EPC_RV', ...
         'all_stn_corr_EPC_RV','all_stn_rmse_EPC_RV','all_stn_CPS','all_stn_CPS_RV', ...
         'var_MRV','var_CPS_RV','var_EPC_RV','var_RVM',...
         'all_stn_corr_CPS','all_stn_rmse_CPS','all_stn_corr_CPS_RV','all_stn_rmse_CPS_RV','-v7.3');
    end

end