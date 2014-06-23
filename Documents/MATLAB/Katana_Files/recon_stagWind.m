% This script will use the proxy groups to reconstruct ENSO using various
% methods, applying staggered windows where possible. It will also plot
% histograms of the methods depending on the plot specified by the variable
% 'theplots'.
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

DIR_NAME = 'Pseudoproxies/glb';  
numstnstocompare = [3 5 10 20 40];                  %%%%%%%%
plotnumber = 1;                                               
for NUM_STNS = numstnstocompare
NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
all_stn_ts=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
all_stn_lat=zeros(NUM_TRIALS,NUM_STNS);
all_stn_lon=zeros(NUM_TRIALS,NUM_STNS);
load([DIR_NAME,'/',num2str(NUM_STNS),'stns_1000prox.mat']);
all_stn_lat=stn_lat;
all_stn_lon=stn_lon;

stn_ts_mn=mean(stn_ts,3); 
stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
for n=1:NUM_TRIALS
    for m=1:NUM_STNS
       all_stn_ts(n,m,:)= (stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m));
    end
end

clear stn_lat stn_lon stn_ts stn_ts_mn stn_ts_std
CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];


%% McGregor et al 2013 method of Median Running Variances

stn_MRV=nan(NUM_TRIALS, NUM_YRS);
stn_corr_MRV = nan(NUM_TRIALS,1);
stn_rmse_MRV = nan(NUM_TRIALS,1);

for n=1:NUM_TRIALS
    stn_movvar=nan(NUM_STNS,NUM_YRS);
    for m=1:NUM_STNS
        stn_movvar(m,:) = movingvar(squeeze(all_stn_ts(n,m,:)),VAR_WDW);
    end
    stn_MRV(n,:) = median(stn_movvar);
end

% Skill Evaluation
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
eof_stn = nan(size(CAL_WDW,1),NUM_TRIALS,NUM_OF_EOFS,NUM_STNS);
PC_stn = nan(size(CAL_WDW,1),NUM_TRIALS,NUM_OF_EOFS,NUM_YRS);
expvar_stn = nan(size(CAL_WDW,1),NUM_TRIALS,NUM_OF_EOFS);
stn_EPC = nan(size(CAL_WDW,1),NUM_TRIALS, NUM_YRS);
stn_EPC_RV = nan(size(CAL_WDW,1),NUM_TRIALS, NUM_YRS);
stn_corr_EPC = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_rmse_EPC = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_corr_EPC_RV = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_rmse_EPC_RV = nan(size(CAL_WDW,1),NUM_TRIALS);

for c=1:size(CAL_WDW,1)
    
for n=1:NUM_TRIALS
    [eof_stn(c,n,:,:),PC_stn(c,n,:,CAL_WDW(c,:)),expvar_stn(c,n,:)] = caleof(squeeze(all_stn_ts(n,:,CAL_WDW(c,:)))', NUM_OF_EOFS, 1);
end

% Flipping EOFs where necessary so PC is similar to PC_stn(1,1,:)
for n=1:NUM_TRIALS
    if corr(n34_ind(CAL_WDW(c,:)), squeeze(PC_stn(c,n,1,CAL_WDW(c,:)))) < 0
        PC_stn(c,n,1,:) = -PC_stn(c,n,1,:);
        eof_stn(c,n,1,:) = -eof_stn(c,n,1,:);
    end
end

% Multiplying EOF and stn data to produce PC timeseries for whole record &
% normalising

for n=1:NUM_TRIALS
    temp = squeeze(eof_stn(c,n,:,:))*squeeze(all_stn_ts(n,:,:));
    stn_EPC(c,n,:) = squeeze(temp(1,:))./std(squeeze(temp(1,:)));
    stn_EPC_RV(c,n,:) = movingvar(squeeze(stn_EPC(c,n,:)),VAR_WDW);
end

% Skill Evaluation

for n=1:NUM_TRIALS
    stn_corr_EPC(c,n) = abs(corr(squeeze(stn_EPC(c,n,:)),n34_ind));
    stn_rmse_EPC(c,n) = sqrt(mean((n34_ind-squeeze(stn_EPC(c,n,:))).^2));
    stn_corr_EPC_RV(c,n) = abs(corr(squeeze(stn_EPC_RV(c,n,RV_WDW)),n34_ind_RV(RV_WDW)));
    stn_rmse_EPC_RV(c,n) = sqrt(mean((n34_ind_RV(RV_WDW)-squeeze(stn_EPC_RV(c,n,RV_WDW))).^2));
end

end
%% Esper et al 2005 CPS Method 

stn_CPS = nan(size(CAL_WDW,1),NUM_TRIALS, NUM_YRS);
stn_CPS_RV = nan(size(CAL_WDW,1),NUM_TRIALS, NUM_YRS);
stn_corr_CPS = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_rmse_CPS = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_corr_CPS_RV = nan(size(CAL_WDW,1),NUM_TRIALS);
stn_rmse_CPS_RV = nan(size(CAL_WDW,1),NUM_TRIALS);

for c=1:size(CAL_WDW,1)

for n=1:NUM_TRIALS
    corr_matrix = corr(n34_ind(CAL_WDW(c,:))*ones(1,NUM_STNS), squeeze(all_stn_ts(n,:,CAL_WDW(c,:)))');
    stn_CPS(c,n,:) = corr_matrix(1,:)*squeeze(all_stn_ts(n,:,:));
end

% Normalising (it already has mean ~0)
for n=1:NUM_TRIALS
    stn_CPS(c,n,:) = squeeze(stn_CPS(c,n,:))./std(squeeze(stn_CPS(c,n,:))');
    stn_CPS_RV(c,n,:) = movingvar(squeeze(stn_CPS(c,n,:)),VAR_WDW);
end

% Skill Evaluation
for n=1:NUM_TRIALS
    stn_corr_CPS(c,n) = corr(squeeze(stn_CPS(c,n,:)),n34_ind);
    stn_rmse_CPS(c,n) = sqrt(mean((n34_ind-squeeze(stn_CPS(c,n,:))).^2));
    stn_corr_CPS_RV(c,n) = corr(squeeze(stn_CPS_RV(c,n,RV_WDW)),n34_ind_RV(RV_WDW));
    stn_rmse_CPS_RV(c,n) = sqrt(mean((n34_ind_RV(RV_WDW)-squeeze(stn_CPS_RV(c,n,RV_WDW))).^2));
end

end
clear c n corr_matrix
%% Plotting

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=28cm y_width=19cm

if strcmp('histline',theplots)
    subplot(length(numstnstocompare),1,plotnumber)
    % Histline Comparisons
    RMSE_XLIM = [0.0 0.5];
    RMSE_YLIM = [0 150];
    CORR_XLIM = [0.0 1.0];
    CORR_YLIM = [0 150];
    CORR_BINS = [CORR_XLIM(1):0.01:CORR_XLIM(2)];
    RMSE_BINS = [RMSE_XLIM(1):0.01:RMSE_XLIM(2)];
    histoline = zeros(size(CORR_BINS,2),size(CAL_WDW,1));
%     histoline = zeros(size(CORR_BINS,2),1); % Use with MRV
    hold on; rgbmap = hsv(size(CAL_WDW,1)); rgbmap(5,3)=0.6;
    for n=1:size(CAL_WDW,1)
        [histoline(:,n) HA] = histline(stn_corr_EPC_RV(n,:),CORR_BINS,'k','plot');    %%%%
%         [histoline HA] = histline(stn_corr_MRV(:),CORR_BINS,'k','plot');  % Use with MRV
        set(HA(1),'Color',rgbmap(n,:))
    end
    hold off; grid on;
    colormap(rgbmap); % colorbar;
    xlim([0 1])
    ylim([0 ceil(max(max(histoline)))])
    title([strrep(DIR_NAME,'_','\_'),' with NUM\_STNS = ',num2str(NUM_STNS)]);
    
%     saveas(gcf,['Plots/recon_stagWind_',DIR_NAME(15:end),'.jpg'])
    if length(numstnstocompare) == plotnumber
        suptitle([strrep(DIR_NAME(15:end),'_','\_'),'\_stn\_corr\_EPC'])
%         saveas(gcf,['Plots/recon_histline_EPC_',DIR_NAME(15:end),'.jpg']);
    end
    plotnumber = plotnumber + 1;
end

end
% Plotting Histogram Comparisons

if strcmp('histogram',theplots)
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
    hist(stn_corr_EPC_RV(1,:),CORR_BINS); % That 1 is the first window
    xlim(CORR_XLIM);
    ylim(CORR_YLIM);
    title('stn\_corr\_EPC\_RV');
    h = findobj(gca,'Type','Patch');
    set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    grid on
    subplot(METH_NUM,2,4)
    hist(stn_rmse_EPC_RV(1,:),RMSE_BINS);
    xlim(RMSE_XLIM);
    ylim(RMSE_YLIM);
    title('stn\_rmse\_EPC\_RV');
    h = findobj(gca,'Type','Patch');
    set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    grid on

    subplot(METH_NUM,2,5) % CPS
    hist(stn_corr_CPS_RV(1,:),CORR_BINS);
    xlim(CORR_XLIM);
    ylim(CORR_YLIM);
    title('stn\_corr\_CPS\_RV');
    h = findobj(gca,'Type','Patch');
    set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    grid on
    subplot(METH_NUM,2,6)
    hist(stn_rmse_CPS_RV(1,:),RMSE_BINS);
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
    hist(stn_rmse_RVM,CORR_BINS);
    xlim(RMSE_XLIM);
    ylim(RMSE_YLIM);
    title('stn\_rmse\_RVM');
    h = findobj(gca,'Type','Patch');
    set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    grid on
    
    suptitle([strrep(DIR_NAME,'_','\_')]);
%     saveas(gcf,['Plots/recon_histogram_',DIR_NAME(15:end),'.jpg'])
end