% This script will use the proxy groups to reconstruct ENSO using a few
% methods. Data needs to be processed from station_select.m
% This script needs mexcdf to be installed in matlab to run, the
% function 'movingvar', 'plotworld', and the folder DataFiles

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';

lat = nc_varget(ts_file,'lat');
lon = nc_varget(ts_file,'lon');
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
%% Loading and normalising proxies

DIR_NAME = 'Pseudoproxies/nplr';
NUM_STNS = 10; NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
all_stn_ts=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
all_stn_lat=zeros(NUM_TRIALS,NUM_STNS);
all_stn_lon=zeros(NUM_TRIALS,NUM_STNS);
for n=1:NUM_TRIALS
    load([DIR_NAME,'/rnd',num2str(NUM_STNS),'_prox',num2str(n),'.mat']);
    all_stn_lat(n,:)=stn_lat;
    all_stn_lon(n,:)=stn_lon;
    stn_ts_mn=mean(stn_ts,2); 
    stn_ts_std=std(stn_ts');
    for m=1:NUM_STNS
	   all_stn_ts(n,m,:)= (stn_ts(m,:)-stn_ts_mn(m))./(stn_ts_std(m));
    end
end
clear stn_lat stn_lon stn_ts stn_ts_mn stn_ts_std

%% McGregor et al 2013 method of Median Running Variances

stn_MRV=nan(NUM_TRIALS, NUM_YRS);
for n=1:NUM_TRIALS
    stn_movvar=nan(NUM_STNS,NUM_YRS);
    for m=1:NUM_STNS
        stn_movvar(m,:) = movingvar(squeeze(all_stn_ts(n,m,:)),VAR_WDW);
    end
    stn_MRV(n,:) = median(stn_movvar);
end

% Skill Evaluation
stn_corr_MRV = nan(NUM_TRIALS,1);
stn_rmse_MRV = nan(NUM_TRIALS,1);
for n=1:NUM_TRIALS
    stn_corr_MRV(n) = corr(squeeze(stn_MRV(n,RV_WDW))',n34_ind_RV(RV_WDW));
    stn_rmse_MRV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_MRV(n,RV_WDW)').^2));
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
    stn_RVM(n,:) = movingvar(med_ts,VAR_WDW);
end

% Skill Evaluation
stn_corr_RVM = nan(NUM_TRIALS,1);
stn_rmse_RVM = nan(NUM_TRIALS,1);
for n=1:NUM_TRIALS
    stn_corr_RVM(n) = corr(squeeze(stn_RVM(n,RV_WDW))',n34_ind_RV(RV_WDW));
    stn_rmse_RVM(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_RVM(n,RV_WDW)').^2));
end
   
%% Braganza et al 2009 method of EOF construction

NUM_OF_EOFS = 3;
eof_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_STNS);
PC_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_YRS);
expvar_stn = nan(NUM_TRIALS,NUM_OF_EOFS);

% EOF calculations using Eigenvector method
for n=1:NUM_TRIALS
    [eof_stn(n,:,:),PC_stn(n,:,:),expvar_stn(n,:)] = caleof(squeeze(all_stn_ts(n,:,:))', NUM_OF_EOFS, 1);
end

% % For checking the correlation values before flipping
% for n=2:NUM_TRIALS
% test_corr(n) = corr(squeeze(PC_stn(1,1,:)), squeeze(PC_stn(n,1,:)));
% end
% plot(sort(test_corr));

% Flipping EOFs where necessary so PC is similar to PC_stn(1,1,:)
for n=2:NUM_TRIALS
    if corr(squeeze(PC_stn(1,1,:)), squeeze(PC_stn(n,1,:))) < 0
        PC_stn(n,1,:) = -PC_stn(n,1,:);
        eof_stn(n,1,:) = -eof_stn(n,1,:);
    end
end

% Normalising (it already has mean 0)
stn_EPC = nan(NUM_TRIALS, NUM_YRS);
stn_EPC_RV = nan(NUM_TRIALS, NUM_YRS);
for n=1:NUM_TRIALS
    stn_EPC(n,:) = squeeze(PC_stn(n,1,:))./std(squeeze(PC_stn(n,1,:))');
    stn_EPC_RV(n,:) = movingvar(stn_EPC(n,:)',VAR_WDW);
end

% Skill Evaluation
stn_corr_EPC = nan(NUM_TRIALS,1);
stn_rmse_EPC = nan(NUM_TRIALS,1);
stn_corr_EPC_RV = nan(NUM_TRIALS,1);
stn_rmse_EPC_RV = nan(NUM_TRIALS,1);
for n=1:NUM_TRIALS
    stn_corr_EPC(n) = abs(corr(squeeze(stn_EPC(n,:))',n34_ind));
    stn_rmse_EPC(n) = sqrt(mean((n34_ind-stn_EPC(n,:)').^2));
    stn_corr_EPC_RV(n) = abs(corr(squeeze(stn_EPC_RV(n,RV_WDW))',n34_ind_RV(RV_WDW)));
    stn_rmse_EPC_RV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_EPC_RV(n,RV_WDW)').^2));
end

%% Esper et al 2005 CPS Method 

CAL_WDW = [300:350];
stn_CPS = nan(NUM_TRIALS, NUM_YRS);
stn_CPS_RV = nan(NUM_TRIALS, NUM_YRS);

for n=1:NUM_TRIALS
    corr_matrix = corr(n34_ind(CAL_WDW)*ones(1,10), squeeze(all_stn_ts(n,:,CAL_WDW))');
    stn_CPS(n,:) = corr_matrix(1,:)*squeeze(all_stn_ts(n,:,:));
end

% Normalising (it already has mean ~0)
for n=1:NUM_TRIALS
    stn_CPS(n,:) = squeeze(stn_CPS(n,:))./std(squeeze(stn_CPS(n,:))');
    stn_CPS_RV(n,:) = movingvar(squeeze(stn_CPS(n,:))',VAR_WDW);
end

% Skill Evaluation
stn_corr_CPS = nan(NUM_TRIALS,1);
stn_rmse_CPS = nan(NUM_TRIALS,1);
stn_corr_CPS_RV = nan(NUM_TRIALS,1);
stn_rmse_CPS_RV = nan(NUM_TRIALS,1);
for n=1:NUM_TRIALS
    stn_corr_CPS(n) = corr(squeeze(stn_CPS(n,:))',n34_ind);
    stn_rmse_CPS(n) = sqrt(mean((n34_ind-stn_CPS(n,:)').^2));
    stn_corr_CPS_RV(n) = corr(squeeze(stn_CPS_RV(n,RV_WDW))',n34_ind_RV(RV_WDW));
    stn_rmse_CPS_RV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_CPS_RV(n,RV_WDW)').^2));
end

%% Plotting Comparisons

RMSE_XLIM = [0.0 0.5];
RMSE_YLIM = [0 150];
CORR_XLIM = [0.0 1.0];
CORR_YLIM = [0 150];
CORR_BINS = [CORR_XLIM(1):0.01:CORR_XLIM(2)];
RMSE_BINS = [RMSE_XLIM(1):0.01:RMSE_XLIM(2)];
METH_NUM = 4;

subplot(METH_NUM,2,1) % MRV
hist(stn_corr_MRV,CORR_BINS);
xlim(CORR_XLIM);
ylim(CORR_YLIM);
title('stn\_corr\_MRV');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on
subplot(METH_NUM,2,2)
hist(stn_rmse_MRV,RMSE_BINS);
xlim(RMSE_XLIM);
ylim(RMSE_YLIM);
title('stn\_rmse\_MRV');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on

subplot(METH_NUM,2,3) % EPC
hist(stn_corr_EPC_RV,CORR_BINS);
xlim(CORR_XLIM);
ylim(CORR_YLIM);
title('stn\_corr\_EPC\_RV');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on
subplot(METH_NUM,2,4)
hist(stn_rmse_EPC_RV,RMSE_BINS);
xlim(RMSE_XLIM);
ylim(RMSE_YLIM);
title('stn\_rmse\_EPC\_RV');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on

subplot(METH_NUM,2,5) % CPS
hist(stn_corr_CPS_RV,CORR_BINS);
xlim(CORR_XLIM);
ylim(CORR_YLIM);
title('stn\_corr\_CPS\_RV');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on
subplot(METH_NUM,2,6)
hist(stn_rmse_CPS_RV,RMSE_BINS);
xlim(RMSE_XLIM);
ylim(RMSE_YLIM);
title('stn\_rmse\_CPS\_RV');
grid on
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');

subplot(METH_NUM,2,7) % RVM
hist(stn_corr_RVM,CORR_BINS);
xlim(CORR_XLIM);
ylim(CORR_YLIM);
title('stn\_corr\_RVM');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on
subplot(METH_NUM,2,8)
hist(stn_rmse_RVM,30);
xlim(RMSE_XLIM);
ylim(RMSE_YLIM);
title('stn\_rmse\_RVM');
h = findobj(gca,'Type','Patch');
set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
grid on