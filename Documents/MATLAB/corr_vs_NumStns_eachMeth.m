% This script will produce a plot of the percentiles of correlation
% corresponding to the number of stations. This is mainly a combination of
% code in recons.m and station_select.m, but modified for many station
% output
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

%% Beggining of Loop

GROUP_NAME = 'glb_ts';
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

load DataFiles/runcorr.mat

NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
lsfrac=nc_varget('DataFiles/sftlf_A1.static.nc','sftlf'); lsfrac(isnan(lsfrac)) = 0 ;

% Selection from areas with absolute correlation over a certain threshold
MIN_COR = 0.3;
STN_MAX = 70;

%% Data Processing and Reconstruction Methods
tic;
STN_NUM_RG = [3:STN_MAX];
all_stn_corr_MRV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_MRV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_corr_RVM = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_RVM = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_EPC = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_corr_EPC = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_EPC_RV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_corr_EPC_RV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_CPS_RV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_corr_CPS_RV = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_rmse_CPS = nan(max(STN_NUM_RG), NUM_TRIALS);
all_stn_corr_CPS = nan(max(STN_NUM_RG), NUM_TRIALS);

for NUM_STNS = STN_NUM_RG
    % Loading and normalising proxies
    load([DIR_NAME,'/CalWdw:1-499/',num2str(NUM_STNS),'stns_1000prox.mat']);
    stn_ts_mn=mean(stn_ts,3);
    stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
    for n=1:NUM_TRIALS
        for m=1:NUM_STNS
           stn_ts(n,m,:)= (stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m));
        end
    end
    
    
    % McGregor et al 2013 method of Median Running Variances
    stn_MRV=nan(NUM_TRIALS, NUM_YRS);
    for n=1:NUM_TRIALS
        stn_movvar=nan(NUM_STNS,NUM_YRS);
        for m=1:NUM_STNS
            stn_movvar(m,:) = movingvar(squeeze(stn_ts(n,m,:)),VAR_WDW);
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
    all_stn_corr_MRV(NUM_STNS,:) = stn_corr_MRV;
    all_stn_rmse_MRV(NUM_STNS,:) = stn_rmse_MRV;
    
    % McGregor et al 2013 method of Running Variance of the Median

    poscorr_all_stn_ts = stn_ts;
    for n=1:NUM_TRIALS
        for m=1:NUM_STNS
            if corr(n34_ind, squeeze(stn_ts(n,m,:))) < 0
                poscorr_all_stn_ts(n,m,:) = -stn_ts(n,m,:);
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
    all_stn_corr_RVM(NUM_STNS,:) = stn_corr_RVM;
    all_stn_rmse_RVM(NUM_STNS,:) = stn_rmse_RVM;
    
    % Braganza et al 2009 method of EOF construction
    NUM_OF_EOFS = 3;
    eof_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_STNS);
    PC_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_YRS);
    expvar_stn = nan(NUM_TRIALS,NUM_OF_EOFS);
    for n=1:NUM_TRIALS
        [eof_stn(n,:,:),PC_stn(n,:,:),expvar_stn(n,:)] = caleof(squeeze(stn_ts(n,:,:))', NUM_OF_EOFS, 1);
    end
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
    all_stn_corr_EPC_RV(NUM_STNS,:) = stn_corr_EPC_RV;
    all_stn_rmse_EPC_RV(NUM_STNS,:) = stn_rmse_EPC_RV;
    all_stn_corr_EPC(NUM_STNS,:) = stn_corr_EPC;
    all_stn_rmse_EPC(NUM_STNS,:) = stn_rmse_EPC;
    
    
    % Esper et al 2005 CPS Method 
    stn_CPS = nan(NUM_TRIALS, NUM_YRS);
    stn_CPS_RV = nan(NUM_TRIALS, NUM_YRS);
    for n=1:NUM_TRIALS
        corr_matrix = corr(n34_ind*ones(1,10), squeeze(stn_ts(n,:,:))');
        stn_CPS(n,:) = corr_matrix(1,:)*squeeze(stn_ts(n,:,:));
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
    all_stn_corr_CPS_RV(NUM_STNS,:) = stn_corr_CPS_RV;
    all_stn_rmse_CPS_RV(NUM_STNS,:) = stn_rmse_CPS_RV;
    all_stn_corr_CPS(NUM_STNS,:) = stn_corr_CPS;
    all_stn_rmse_CPS(NUM_STNS,:) = stn_rmse_CPS;
    toc;
end % Took 40 min for 70 

save('DataFiles/500yrCalWdw_meth_stats.mat','all_stn_corr_MRV','all_stn_rmse_MRV', ...
     'all_stn_corr_RVM','all_stn_corr_RVM', ...
     'all_stn_corr_EPC', 'all_stn_rmse_EPC', ...
     'all_stn_corr_EPC_RV', 'all_stn_rmse_EPC_RV', ...
     'all_stn_corr_CPS', 'all_stn_rmse_CPS', ...
     'all_stn_corr_CPS_RV', 'all_stn_rmse_CPS_RV'         );

% %% Making the plot
% 
% qEPC = quantile(all_stn_corr_EPC_RV,[.05 .5 .95], 2);
% qMRV = quantile(all_stn_corr_MRV,[.05 .5 .95], 2);
% qCPS = quantile(all_stn_corr_CPS_RV,[.05 .5 .95], 2);
% 
% clf; axes; hold on; grid on;
% Hnd(1,1) = plot(squeeze(qMRV(:,1)),'--g');
% set(Hnd(1,1),'Color','g','MarkerFaceColor','g');
% Hnd(1,3) = plot(squeeze(qMRV(:,3)),'--g');
% Hnd(2,1) = plot(squeeze(qEPC(:,1)),'--r');
% Hnd(2,3) = plot(squeeze(qEPC(:,3)),'--r');
% Hnd(3,1) = plot(squeeze(qCPS(:,1)),'--b');
% Hnd(3,3) = plot(squeeze(qCPS(:,3)),'--b');
% Hnd(1,2) = plot(squeeze(qMRV(:,2)),'g','LineWidth',2);
% Hnd(2,2) = plot(squeeze(qEPC(:,2)),'r','LineWidth',2);
% Hnd(3,2) = plot(squeeze(qCPS(:,2)),'b','LineWidth',2);
% hold off;
% set(Hnd(1,[1,3]),'Color','g','MarkerFaceColor','g');
% set(Hnd(2,[1,3]),'Color','r','MarkerFaceColor','r');
% set(Hnd(3,[1,3]),'Color','b','MarkerFaceColor','b');
% ylabel('Percentile Correlations with Nino3.4 index','FontSize',14  ); 
% xlabel('Number of Stations used in reconstruction','FontSize',14  );
% legend([Hnd(1:3,2); Hnd(3,1); Hnd(3,3);],'MRV Median','EPC\_RV Median','CPS\_RV Median', ...
%        '5^t^h Percentile','95^t^h Percentile','location','southeast'                );
% % title('Reconstructions from a global selection of pseudoproxies, using all 499 years of data')
% set(gca, ...
%     'TickDir','in', ...
%     'Box', 'on',    ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'on'      , ...
%     'LineWidth'   , 1,  ...    
%     'YLim'        , [0 1], ...
%     'XLim'        , [0 70], ...
%     'FontSize'    , 14 ...
%     );
% 
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm
% saveas(gcf,['Plots/corr_vsNumStns_eachMeth_glb.jpg'])