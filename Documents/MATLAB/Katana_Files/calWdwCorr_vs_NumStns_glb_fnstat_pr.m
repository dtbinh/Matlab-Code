% This script will calculate the range of the error in the calibration windows 
% for varying station numbers, using staggered windows. It also outputs a
% mat file containing alot of statistics produced for each method,
% num_stns, and calibration window. Series 2

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
% window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];

%% Loading proxies

DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/glb_pr_fnstat'];                                 %%%%%%%
numstnstocompare = [3:70];                                     %%%%%%%%
NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
% Calibration windows set to being 10 overlapping windows over 499 years
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
% CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];

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

all_stn_pr=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat']);
all_stn_lat=stn_lat;
all_stn_lon=stn_lon;

stn_pr_mn=mean(stn_pr,3); 
stn_pr_std=squeeze(std(permute(stn_pr,[3,1,2])));
for n=1:NUM_TRIALS
    for m=1:NUM_STNS
       all_stn_pr(n,m,:)= single((stn_pr(n,m,:)-stn_pr_mn(n,m))./(stn_pr_std(n,m)));
    end
end

clear stn_lat stn_lon stn_pr stn_pr_mn stn_pr_std

%% McGregor et al 2013 method of Median Running Variances

stn_MRV=nan(NUM_TRIALS, NUM_YRS);
stn_corr_MRV = nan(NUM_TRIALS,1);
stn_rmse_MRV = nan(NUM_TRIALS,1);

for n=1:NUM_TRIALS
    stn_movvar=nan(NUM_STNS,NUM_YRS);
    for m=1:NUM_STNS
        stn_movvar(m,:) = movingvar(squeeze(all_stn_pr(n,m,:)),VAR_WDW);
    end
    stn_MRV(n,:) = single(median(stn_movvar));
end

% Skill Evaluation
for n=1:NUM_TRIALS
    stn_corr_MRV(n) = single(corr(squeeze(stn_MRV(n,RV_WDW))',n34_ind_RV(RV_WDW)));
    stn_rmse_MRV(n) = single(sqrt(mean((n34_ind_RV(RV_WDW)-stn_MRV(n,RV_WDW)').^2)));
end

%% McGregor et al 2013 method of Running Variance of the Median

poscorr_all_stn_pr = all_stn_pr;
for n=1:NUM_TRIALS
    for m=1:NUM_STNS
        if corr(n34_ind, squeeze(all_stn_pr(n,m,:))) < 0
            poscorr_all_stn_pr(n,m,:) = -all_stn_pr(n,m,:);
        end
    end
end

stn_RVM=nan(NUM_TRIALS, NUM_YRS);
for n=1:NUM_TRIALS
    med_pr = median(squeeze(poscorr_all_stn_pr(n,:,:)))';
    med_pr = (med_pr-mean(med_pr))./std(med_pr);
    stn_RVM(n,:) = single(movingvar(med_pr,VAR_WDW));
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
    [eof_stn(n,:,:),PC_stn(n,:,CAL_WDW(c,:)),expvar_stn(n,:)] = caleof(squeeze(all_stn_pr(n,:,CAL_WDW(c,:)))', NUM_OF_EOFS, 1);
end

% Flipping EOFs where necessary so PC is similar to Nino34 correlation
for n=1:NUM_TRIALS
    if corr(n34_ind(CAL_WDW(c,:)), squeeze(PC_stn(n,1,CAL_WDW(c,:)))) < 0
        PC_stn(n,1,:) = -PC_stn(n,1,:);
        eof_stn(n,1,:) = -eof_stn(n,1,:);
    end
end

% Multiplying EOF and stn data to produce PC timeseries for whole record &
% normalising

for n=1:NUM_TRIALS
    temp = squeeze(eof_stn(n,:,:))*squeeze(all_stn_pr(n,:,:));
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
    corr_matrix = corr(n34_ind(CAL_WDW(c,:))*ones(1,NUM_STNS), squeeze(all_stn_pr(n,:,CAL_WDW(c,:)))');
    stn_CPS(n,:) = corr_matrix(1,:)*squeeze(all_stn_pr(n,:,:));
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
     'all_stn_corr_CPS','all_stn_rmse_CPS','all_stn_corr_CPS_RV','all_stn_rmse_CPS_RV','-v7.3','window');

end

% %% Plotting Method
% DIR_NAME = 'Pseudoproxies/glb';
% CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];
% temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
% temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
% temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
% 
% for c=1:size(CAL_WDW,1)
%     load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
%      'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV')
%     temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
%     temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
%     temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
% end
% method = 'all';
% 
% % Plotting EPC
% 
% if strcmp(method,'EPC') | strcmp(method,'all')
% subplot(1,3,1)
%     corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% %     clf; axes; hold on; rgbmap = jet(size(CAL_WDW,1)); rgbmap(2:3,:) = rgbmap(1:2,:); rgbmap(1,:) = [0 0 0];
% %     HA = nan(size(CAL_WDW,1),3);
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,1) = plot(squeeze(corr_RV_qn(:,c,1)),'-v');
% %         set(HA(c,1),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %         HA(c,2) = plot(squeeze(corr_RV_qn(:,c,3)),'-^');
% %         set(HA(c,2),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %     end
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,3) = plot(squeeze(corr_RV_qn(:,c,2)),'LineWidth',2);
% %         set(HA(c,3),'Color',rgbmap(c,:));
% %     end
% %     title([strrep(DIR_NAME(15:end),'_','\_'), ' - Correlation percentiles of EPC\_RV for each calibration window'])
% %     xlabel('Number of Stations included in reconstruction'); ylabel('Correlation')
% %     grid on
% %     ylim([0 1])
% %     legend([HA(1:end,3);HA(1,1);HA(1,2)],'Year 1-50 median','Year 51-100 median','Year 101-150 median','Year 151-200 median',...
% %                        'Year 201-250 median','Year 251-300 median','Year 301-350 median','Year 351-400 median',...
% %                        'Year 401-450 median','Year 450-499 median', '5^t^h Percentile','95^t^h Percentile',...
% %                        'location','southeast' );
% 
% % Range Plotting
%     
%     corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
%     corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
%     corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','g','k','add',0.5);
% %     legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast')
%     xlim([0,70]); ylim([0,1]); grid on
%     ylabel('Correlation')
%     title(['EPC\_RV'])
%     
% %     saveas(gcf,['Plots/calWdwCorr_vs_NumStns_',DIR_NAME(15:end),'_epc.jpg'])
% end
% 
% if strcmp(method,'CPS') | strcmp(method,'all') 
% % Plotting CPS
% subplot(1,3,2)
%     corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
% %     clf; axes; hold on; rgbmap = jet(size(CAL_WDW,1)); rgbmap(2:3,:) = rgbmap(1:2,:); rgbmap(1,:) = [0 0 0];
% %     HA = nan(size(CAL_WDW,1),3);
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,1) = plot(squeeze(corr_RV_qn(:,c,1)),'-v');
% %         set(HA(c,1),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %         HA(c,2) = plot(squeeze(corr_RV_qn(:,c,3)),'-^');
% %         set(HA(c,2),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %     end
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,3) = plot(squeeze(corr_RV_qn(:,c,2)),'LineWidth',2);
% %         set(HA(c,3),'Color',rgbmap(c,:));
% %     end
% %     xlabel('Number of Stations included in reconstruction'); ylabel('Correlation')
% %     grid on
% %     ylim([0 1])
% %     legend([HA(1:end,3);HA(1,1);HA(1,2)],'Year 1-50 median','Year 51-100 median','Year 101-150 median','Year 151-200 median',...
% %                        'Year 201-250 median','Year 251-300 median','Year 301-350 median','Year 351-400 median',...
% %                        'Year 401-450 median','Year 450-499 median', '5^t^h Percentile','95^t^h Percentile',...
% %                        'location','southeast' );
%                    
%     corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
%     corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
%     corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','g','k','add',0.5);
% %     legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast')
%     xlim([0,70]); ylim([0,1]); grid on
%     xlabel('Number of Stations included in reconstruction');
%     
%     title(['CPS\_RV'])
% %     saveas(gcf,['Plots/calWdwCorr_vs_NumStns_',DIR_NAME(15:end),'_cps.jpg'])
% end
%     
% if strcmp(method,'MRV') | strcmp(method,'all')
% % Plotting MRV
% subplot(1,3,3)
%     corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
% %     clf; axes; hold on; rgbmap = jet(size(CAL_WDW,1)); rgbmap(2:3,:) = rgbmap(1:2,:); rgbmap(1,:) = [0 0 0];
% %     HA = nan(size(CAL_WDW,1),3);
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,1) = plot(squeeze(corr_RV_qn(:,c,1)),'-v');
% %         set(HA(c,1),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %         HA(c,2) = plot(squeeze(corr_RV_qn(:,c,3)),'-^');
% %         set(HA(c,2),'Color',rgbmap(c,:),'MarkerFaceColor',rgbmap(c,:));
% %     end
% %     for c=1:size(CAL_WDW,1)
% %         HA(c,3) = plot(squeeze(corr_RV_qn(:,c,2)),'LineWidth',2);
% %         set(HA(c,3),'Color',rgbmap(c,:));
% %     end
% %     title([strrep(DIR_NAME(15:end),'_','\_'), ' - Correlation percentiles of MRV for each calibration window proxy group'])
% %     xlabel('Number of Stations included in reconstruction'); ylabel('Correlation')
% %     grid on
% %     ylim([0 1])
% %     legend([HA(1:end,3);HA(1,1);HA(1,2)],'Year 1-50 median','Year 51-100 median','Year 101-150 median','Year 151-200 median',...
% %                        'Year 201-250 median','Year 251-300 median','Year 301-350 median','Year 351-400 median',...
% %                        'Year 401-450 median','Year 450-499 median', '5^t^h Percentile','95^t^h Percentile',...
% %                        'location','southeast' );
% %     clf;
% %     HA(1) = plot(squeeze(corr_RV_qn(:,1)),'-kv'); hold on;
% %     set(HA(1),'MarkerFaceColor','k');
% %     HA(2) = plot(squeeze(corr_RV_qn(:,3)),'-k^');
% %     set(HA(2),'MarkerFaceColor','k');
% %     HA(3) = plot(squeeze(corr_RV_qn(:,2)),'k','LineWidth',2); hold off;
% %
% %     xlabel('Number of Stations included in reconstruction'); ylabel('Correlation')
% %     grid on
% %     ylim([0 1])
% %     legend('5^t^h Percentile','95^t^h Percentile','Median',...
% %                        'location','southeast' );
% 
%     corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
%     corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
%     corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
%     jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','g','k','add',0.5);
%     legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast')
%     xlim([0,70]); ylim([0,1]); grid on
% %     xlabel('Number of Stations included in reconstruction'); ylabel('Correlation')
%         
%     title(['MRV'])
% end
% 
% if strcmp(method,'all')
% %     close
%     suptitle([strrep(DIR_NAME(15:end),'_','\_'),' - Ranges of Correlation percentiles'])
%     legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast');
%     set(gcf, 'PaperPosition', [0 0 28 19]);
%     for i=1:3
% %         subplot(1,3,i)
%         set(gca, 'FontSize',14, ...
%             'LineWidth', 1.0, ...
%             'Box', 'on', ...
%             'YTick', [0:0.1:1]); 
%     end
%         
%         
%         
%     set(legendH, 'FontSize',10);
%     saveas(gcf,['Plots/calWdwCorr_vs_NumStns_',DIR_NAME(15:end),'.jpg'])
% end