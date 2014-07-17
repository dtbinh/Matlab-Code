% This script will calculate the range of the error in the calibration windows 
% for varying station numbers, using staggered windows. It also outputs a
% mat file containing alot of statistics produced for each method,
% num_stns, and calibration window.

% Should look like 10 percentile lines for each calibration window, vs the num_stns.

% This script needs mexcdf to be installed in matlab to run, the
% function 'movingvar', 'plotworld', and the folder DataFiles

%% Setup
tic;
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

%% Loading proxies
GROUP_NAME = 'glb_ts';
% DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
numstnstocompare = [3:70];                                     %%%%%%%%
NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
% Calibration windows set to being 10 overlapping windows over 499 years
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

%% Beginning the Loop
for c=1:size(CAL_WDW,1)

all_stn_MRV=nan(max(numstnstocompare), NUM_TRIALS, NUM_YRS,'single');
all_stn_corr_MRV = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_rmse_MRV = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_EPC = nan(max(numstnstocompare),NUM_TRIALS, NUM_YRS,'single');
all_stn_EPC_RV = nan(max(numstnstocompare),NUM_TRIALS, NUM_YRS,'single');
all_stn_corr_EPC = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_rmse_EPC = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_corr_EPC_RV = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_rmse_EPC_RV = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_CPS = nan(max(numstnstocompare),NUM_TRIALS, NUM_YRS,'single');
all_stn_CPS_RV = nan(max(numstnstocompare),NUM_TRIALS, NUM_YRS,'single');
all_stn_corr_CPS = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_rmse_CPS = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_corr_CPS_RV = nan(max(numstnstocompare),NUM_TRIALS,'single');
all_stn_rmse_CPS_RV = nan(max(numstnstocompare),NUM_TRIALS,'single');
    
    
for NUM_STNS = numstnstocompare

all_stn_ts=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat']);
all_stn_lat=stn_lat;
all_stn_lon=stn_lon;

stn_ts_mn=mean(stn_ts,3); 
stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
for n=1:NUM_TRIALS
    for m=1:NUM_STNS
       all_stn_ts(n,m,:)= single((stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m)));
    end
end

clear stn_lat stn_lon stn_ts stn_ts_mn stn_ts_std

%% McGregor et al 2013 method of Median Running Variances

stn_MRV=nan(NUM_TRIALS, NUM_YRS);
stn_corr_MRV = nan(NUM_TRIALS,1);
stn_rmse_MRV = nan(NUM_TRIALS,1);

for n=1:NUM_TRIALS
    stn_movvar=nan(NUM_STNS,NUM_YRS);
    for m=1:NUM_STNS
        stn_movvar(m,:) = movingvar(squeeze(all_stn_ts(n,m,:)),VAR_WDW);
    end
    stn_MRV(n,:) = single(median(stn_movvar));
end

% Skill Evaluation
for n=1:NUM_TRIALS
    stn_corr_MRV(n) = single(corr(squeeze(stn_MRV(n,RV_WDW))',n34_ind_RV(RV_WDW)));
    stn_rmse_MRV(n) = single(sqrt(mean((n34_ind_RV(RV_WDW)-stn_MRV(n,RV_WDW)').^2)));
end

%% McGregor et al 2013 method of Running Variance of the Median

poscorr_all_stn_ts = all_stn_ts;
for n=1:NUM_TRIALS
    for m=1:NUM_STNS
        if corr(n34_ind, squeeze(all_stn_ts(n,m,:))) < 0
            poscorr_all_stn_ts(n,m,:) = -all_stn_ts(n,m,:);
        end
    end
end

stn_RVM=nan(NUM_TRIALS, NUM_YRS);
for n=1:NUM_TRIALS
    med_ts = median(squeeze(poscorr_all_stn_ts(n,:,:)))';
    med_ts = (med_ts-mean(med_ts))./std(med_ts);
    stn_RVM(n,:) = single(movingvar(med_ts,VAR_WDW));
end

% Skill Evaluation
stn_corr_RVM = nan(NUM_TRIALS,1);
stn_rmse_RVM = nan(NUM_TRIALS,1);
for n=1:NUM_TRIALS
    stn_corr_RVM(n) = single(corr(squeeze(stn_RVM(n,RV_WDW))',n34_ind_RV(RV_WDW)));
    stn_rmse_RVM(n) = single(sqrt(mean((n34_ind_RV(RV_WDW)-stn_RVM(n,RV_WDW)').^2)));
end
%% Braganza et al 2009 method of EOF construction

NUM_OF_EOFS = 3;
eof_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_STNS);
PC_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_YRS);
expvar_stn = nan(NUM_TRIALS,NUM_OF_EOFS);
stn_EPC = nan(NUM_TRIALS, NUM_YRS);
stn_EPC_RV = nan(NUM_TRIALS, NUM_YRS);
stn_corr_EPC = nan(NUM_TRIALS,1);
stn_rmse_EPC = nan(NUM_TRIALS,1);
stn_corr_EPC_RV = nan(NUM_TRIALS,1);
stn_rmse_EPC_RV = nan(NUM_TRIALS,1);

for n=1:NUM_TRIALS
    [eof_stn(n,:,:),PC_stn(n,:,CAL_WDW(c,:)),expvar_stn(n,:)] = caleof(squeeze(all_stn_ts(n,:,CAL_WDW(c,:)))', NUM_OF_EOFS, 1);
end

% Flipping EOFs where necessary so PC is similar to PC_stn(1,1,:)
for n=1:NUM_TRIALS
    if corr(n34_ind(CAL_WDW(c,:)), squeeze(PC_stn(n,1,CAL_WDW(c,:)))) < 0
        PC_stn(n,1,:) = -PC_stn(n,1,:);
        eof_stn(n,1,:) = -eof_stn(n,1,:);
    end
end

% Multiplying EOF and stn data to produce PC timeseries for whole record &
% normalising

for n=1:NUM_TRIALS
    temp = squeeze(eof_stn(n,:,:))*squeeze(all_stn_ts(n,:,:));
    stn_EPC(n,:) = single(squeeze(temp(1,:))./std(squeeze(temp(1,:))));
    stn_EPC_RV(n,:) = single(movingvar(squeeze(stn_EPC(n,:)'),VAR_WDW));
end

% Skill Evaluation

for n=1:NUM_TRIALS
    stn_corr_EPC(n) = single(abs(corr(squeeze(stn_EPC(n,:))',n34_ind)));
    stn_rmse_EPC(n) = single(sqrt(mean((n34_ind-squeeze(stn_EPC(n,:))').^2)));
    stn_corr_EPC_RV(n) = single(abs(corr(squeeze(stn_EPC_RV(n,RV_WDW))',n34_ind_RV(RV_WDW))));
    stn_rmse_EPC_RV(n) = single(sqrt(mean((n34_ind_RV(RV_WDW)-squeeze(stn_EPC_RV(n,RV_WDW))').^2)));
end


%% Esper et al 2005 CPS Method

stn_CPS = nan(NUM_TRIALS, NUM_YRS);
stn_CPS_RV = nan(NUM_TRIALS, NUM_YRS);
stn_corr_CPS = nan(NUM_TRIALS,1);
stn_rmse_CPS = nan(NUM_TRIALS,1);
stn_corr_CPS_RV = nan(NUM_TRIALS,1);
stn_rmse_CPS_RV = nan(NUM_TRIALS,1);

for n=1:NUM_TRIALS
    corr_matrix = corr(n34_ind(CAL_WDW(c,:))*ones(1,NUM_STNS), squeeze(all_stn_ts(n,:,CAL_WDW(c,:)))');
    stn_CPS(n,:) = corr_matrix(1,:)*squeeze(all_stn_ts(n,:,:));
end

% Normalising (it already has mean ~0)
for n=1:NUM_TRIALS
    stn_CPS(n,:) = single(squeeze(stn_CPS(n,:))./std(squeeze(stn_CPS(n,:))'));
    stn_CPS_RV(n,:) = single(movingvar(squeeze(stn_CPS(n,:))',VAR_WDW));
end

% Skill Evaluation
for n=1:NUM_TRIALS
    stn_corr_CPS(n) = single(corr(squeeze(stn_CPS(n,:))',n34_ind));
    stn_rmse_CPS(n) = single(sqrt(mean((n34_ind-squeeze(stn_CPS(n,:))').^2)));
    stn_corr_CPS_RV(n) = single(corr(squeeze(stn_CPS_RV(n,RV_WDW))',n34_ind_RV(RV_WDW)));
    stn_rmse_CPS_RV(n) = single(sqrt(mean((n34_ind_RV(RV_WDW)-squeeze(stn_CPS_RV(n,RV_WDW))').^2)));
end

all_stn_MRV(NUM_STNS,:,:) = stn_MRV;
all_stn_corr_MRV(NUM_STNS,:) = stn_corr_MRV;
all_stn_rmse_MRV(NUM_STNS,:) = stn_rmse_MRV;
all_stn_EPC(NUM_STNS,:,:) = stn_EPC;
all_stn_EPC_RV(NUM_STNS,:,:) = stn_EPC_RV;
all_stn_corr_EPC(NUM_STNS,:) = stn_corr_EPC;
all_stn_rmse_EPC(NUM_STNS,:) = stn_rmse_EPC;
all_stn_corr_EPC_RV(NUM_STNS,:) = stn_corr_EPC_RV;
all_stn_rmse_EPC_RV(NUM_STNS,:) = stn_rmse_EPC_RV;
all_stn_CPS(NUM_STNS,:,:) = stn_CPS;
all_stn_CPS_RV(NUM_STNS,:,:) = stn_CPS_RV;
all_stn_corr_CPS(NUM_STNS,:) = stn_corr_CPS;
all_stn_rmse_CPS(NUM_STNS,:) = stn_rmse_CPS;
all_stn_corr_CPS_RV(NUM_STNS,:) = stn_corr_CPS_RV;
all_stn_rmse_CPS_RV(NUM_STNS,:) = stn_rmse_CPS_RV;
toc;
end
% That took about an hour per calibration window
save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
     'all_stn_MRV','all_stn_corr_MRV','all_stn_rmse_MRV', ...
     'DIR_NAME','CAL_WDW','all_stn_EPC','all_stn_corr_EPC','all_stn_rmse_EPC','all_stn_EPC_RV', ...
     'all_stn_corr_EPC_RV','all_stn_rmse_EPC_RV','all_stn_CPS','all_stn_CPS_RV', ...
     'all_stn_corr_CPS','all_stn_rmse_CPS','all_stn_corr_CPS_RV','all_stn_rmse_CPS_RV','-v7.3');

end
