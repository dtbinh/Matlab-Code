% This script will determine statistical significance of a reconstruction's correlation to Nino3.4 Index

%% Setup
tic;
load DataFiles/model_output.mat

VAR_WDW = 30; % Moving window for moving variance is 30 Years
window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)]; numstnstocompare = 3:70; NUM_TRIALS=1000; NUM_YRS=499;

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
% DIR_NAME = ['/home/nfs/z3372730/Documents/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

temp_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,NUM_YRS,'single');
temp_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,NUM_YRS,'single');
temp_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,NUM_YRS,'single');
temp_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,NUM_YRS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_EPC_RV','all_stn_CPS_RV','all_stn_MRV','all_stn_RVM')
    temp_EPC_RV(:,c,:,:) = all_stn_EPC_RV;
    temp_CPS_RV(:,c,:,:) = all_stn_CPS_RV;
    temp_MRV(:,c,:,:) = all_stn_MRV;
    temp_RVM(:,c,:,:) = all_stn_RVM;
end

all_corr = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_sig_level = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_df = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
for n=numstnstocompare
    for c=1:size(CAL_WDW,1)
        for i=1:NUM_TRIALS
            one_recon = squeeze(temp_MRV(n,c,i,:));
            one_recon = one_recon - mean(one_recon);
            [corr_coef,sig_level,degrees_of_freedom] = calc_statsig(n34_ind_RV(RV_WDW),one_recon(RV_WDW));
            all_corr(n,c,i) = corr_coef;
            all_corr_sig_level(n,c,i) = sig_level;
            all_corr_df(n,c,i) = degrees_of_freedom;
        end
    end
end

save('../../Dropbox/sig_corrs_glb_ts_MRV.mat','all_corr','all_corr_sig_level','all_corr_df');


all_corr = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_sig_level = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_df = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
for n=numstnstocompare
    for c=1:size(CAL_WDW,1)
        for i=1:NUM_TRIALS
            one_recon = squeeze(temp_RVM(n,c,i,:));
            one_recon = one_recon - mean(one_recon);
            [corr_coef,sig_level,degrees_of_freedom] = calc_statsig(n34_ind_RV(RV_WDW),one_recon(RV_WDW));
            all_corr(n,c,i) = corr_coef;
            all_corr_sig_level(n,c,i) = sig_level;
            all_corr_df(n,c,i) = degrees_of_freedom;
        end
    end
end

save('../../Dropbox/sig_corrs_glb_ts_RVM.mat','all_corr','all_corr_sig_level','all_corr_df');


all_corr = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_sig_level = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_df = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
for n=numstnstocompare
    for c=1:size(CAL_WDW,1)
        for i=1:NUM_TRIALS
            one_recon = squeeze(temp_CPS_RV(n,c,i,:));
            one_recon = one_recon - mean(one_recon);
            [corr_coef,sig_level,degrees_of_freedom] = calc_statsig(n34_ind_RV(RV_WDW),one_recon(RV_WDW));
            all_corr(n,c,i) = corr_coef;
            all_corr_sig_level(n,c,i) = sig_level;
            all_corr_df(n,c,i) = degrees_of_freedom;
        end
    end
end

save('../../Dropbox/sig_corrs_glb_ts_CPS_RV.mat','all_corr','all_corr_sig_level','all_corr_df');


all_corr = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_sig_level = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
all_corr_df = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
for n=numstnstocompare
    for c=1:size(CAL_WDW,1)
        for i=1:NUM_TRIALS
            one_recon = squeeze(temp_EPC_RV(n,c,i,:));
            one_recon = one_recon - mean(one_recon);
            [corr_coef,sig_level,degrees_of_freedom] = calc_statsig(n34_ind_RV(RV_WDW),one_recon(RV_WDW));
            all_corr(n,c,i) = corr_coef;
            all_corr_sig_level(n,c,i) = sig_level;
            all_corr_df(n,c,i) = degrees_of_freedom;
        end
    end
end

save('../../Dropbox/sig_corrs_glb_ts_EPC_RV.mat','all_corr','all_corr_sig_level','all_corr_df');


toc;
%% Autocorrelation of GFDL CM2.1 Data
% 
% 
% [acf, lags,bounds] = autocorr(n34_ind_RV-mean(n34_ind_RV),200);
% plot(lags,acf); xlabel('Lag (yrs)'); ylabel('Auto-Correlation');
% % 

