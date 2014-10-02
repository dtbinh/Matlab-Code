%% This script will compare the RVM and MRV methods in the reconstructions

%% Setup
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
RV_WDW = [15:(499-14)]; NUM_YRS=499; NUM_TRIALS=1000; NUM_STNS=10;
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
numstnstocompare=3:70;

% %% Loading Options
% c=1;
% GROUP_NAME = 'glb_ts';
% DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
% load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox'])
% reconstruction_iteration = 104;
% %% Calculations
% % Normalising
% stn_ts_nm=nan(size(stn_ts)); % nm means normalised
% stn_ts_mn=mean(stn_ts,3); 
% stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
% for n=1:NUM_TRIALS
%     for m=1:NUM_STNS
%        stn_ts_nm(n,m,:)= single((stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m)));
%     end
% end
% 
% poscorr_stn_ts = stn_ts_nm;
% for n=1:NUM_TRIALS
%     for m=1:NUM_STNS
%         if corr(n34_ind, squeeze(stn_ts_nm(n,m,:))) < 0
%             poscorr_stn_ts(n,m,:) = -stn_ts_nm(n,m,:);
%         end
%     end
% end
% % Calculating Running variance
% stn_movvar=nan(NUM_TRIALS,NUM_STNS,NUM_YRS);
% for n=1:NUM_TRIALS
%     for m=1:NUM_STNS
%         stn_movvar(n,m,:) = movingvar(squeeze(stn_ts_nm(n,m,:)),VAR_WDW);
%     end
% end
% stn_MRV = squeeze(single(median(permute(stn_movvar,[2,1,3]))));
% 
% stn_RVM = nan(size(stn_MRV));
% med_ts_nm = nan(size(stn_MRV));
% for n=1:NUM_TRIALS
%     med_ts = median(squeeze(poscorr_stn_ts(n,:,:)))';
%     med_ts_nm(n,:) = (med_ts-mean(med_ts))./std(med_ts)';
%     stn_RVM(n,:) = single(movingvar(squeeze(med_ts_nm(n,:)'),VAR_WDW));
% end
% 
% %% Plotting
% mov_corr = nan(NUM_STNS,NUM_YRS);
% for n=1:5
%     subplot(5,1,n)
% %     plot(n34_ind,'k','LineWidth',2); % 104 looks bad
%     mov_corr(n,:) = movingCorrelation([n34_ind,squeeze(poscorr_stn_ts(reconstruction_iteration,n,:))],VAR_WDW,2);
%     area(squeeze(mov_corr(n,32:end)));
%     xlim([0 500]); ylim([-1 1]); grid on
%     title(['Lon:',num2str(lon(stn_lon(reconstruction_iteration,n))),'Lat:',num2str(lat(stn_lat(reconstruction_iteration,n)))])
% %     subplot(5,2,2*n)
% %     plot(squeeze(stn_movvar(104,n,:)),'b') % 104 looks bad
% %     hold on; plot(n34_ind_RV,'k'); hold off;
% %     xlim([0 500]); ylim([0 3]); grid on
%     
% end
% figure
% for n=6:10
%     subplot(5,1,n-5)
% %     plot(n34_ind,'k','LineWidth',2); % 104 looks bad
%     mov_corr(n,:) = movingCorrelation([n34_ind,squeeze(poscorr_stn_ts(reconstruction_iteration,n,:))],VAR_WDW,2);
%     area(squeeze(mov_corr(n,32:end))); xlim([0 500]); ylim([-1 1]); grid on
%     title(['Lon:',num2str(lon(stn_lon(reconstruction_iteration,n))),'Lat:',num2str(lat(stn_lat(reconstruction_iteration,n)))])
%    
% %     subplot(5,2,2*n)
% %     plot(squeeze(stn_movvar(104,n,:)),'b') % 104 looks bad
% %     hold on; plot(n34_ind_RV,'k'); hold off;
% %     xlim([0 500]); ylim([0 3]); grid on
%     
% end
% 
% figure
% % Running Variance
% for n=1:5
%     subplot(5,1,n)
% %     plot(n34_ind,'k','LineWidth',2); % 104 looks bad
% %     hold on; plot(squeeze(stn_ts_nm(reconstruction_iteration,n,:)),'r'); hold off;
% %     xlim([0 500]); ylim([-3 3]); grid on
% %     subplot(5,2,2*n)
%     plot(squeeze(stn_movvar(104,n,:)),'b') % 104 looks bad
%     hold on; plot(n34_ind_RV,'k'); hold off;
%     xlim([0 500]); ylim([0 3]); grid on
%     title(['Lon:',num2str(lon(stn_lon(reconstruction_iteration,n))),'Lat:',num2str(lat(stn_lat(reconstruction_iteration,n)))])
%  
%     
% end
% figure
% for n=6:10
%     subplot(5,1,n-5)
% %     plot(n34_ind,'k','LineWidth',2); % 104 looks bad
% %     hold on; plot(squeeze(stn_ts_nm(reconstruction_iteration,n,:)),'r'); hold off;
% %     xlim([0 500]); ylim([-3 3]); grid on
% %     subplot(5,2,2*n)
%     plot(squeeze(stn_movvar(104,n,:)),'b') % 104 looks bad
%     hold on; plot(n34_ind_RV,'k'); hold off;
%     xlim([0 500]); ylim([0 3]); grid on
%     title(['Lon:',num2str(lon(stn_lon(reconstruction_iteration,n))),'Lat:',num2str(lat(stn_lat(reconstruction_iteration,n)))])
% 
% end
% 
% figure
% % Plotting Medians etc
% subplot(4,1,1)
% plot(n34_ind,'k','LineWidth',2); hold on;
% plot(squeeze(med_ts_nm(reconstruction_iteration,:)),'r'); hold off;
% title('Normalised Median of Time series')
% 
% subplot(4,1,2)
% plot(n34_ind_RV,'k','LineWidth',2); hold on;
% plot(squeeze(stn_RVM(reconstruction_iteration,:)),'r'); hold off;
% title('RVM method')
% 
% subplot(4,1,3)
% plot(n34_ind_RV,'k','LineWidth',2); hold on;
% plot(squeeze(stn_MRV(reconstruction_iteration,:)),'r'); hold off;
% title('MRV method')
% 
% subplot(4,1,4)
% temp=std(mov_corr); temp=temp(31:end);
% rvm=squeeze(stn_RVM(reconstruction_iteration,:)); rvm=rvm(16:end-15);
% [axH std_H rvm_H] = plotyy(1:469,temp,1:469,rvm);
% legend([std_H rvm_H],'std(r(ts))','RVM','orientation','horizontal')
% xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
% title(['std of corr(proxy temp) vs RVM method, r=',num2str(corr(temp',rvm'))])
% 
% % Scatterplot
% 
% scatter(rvm,temp)
% xlabel('rvm'); ylabel('std(r(ts))');


%% Scatterplot - Reconstruction skill vs correlation to std(corr(stn_ts,n34_ind))
tic;
window = 31; % The running window in years
GROUP_NAME = 'glb_ts';

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

EPC_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
CPS_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
MRV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
RVM_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
EPC_RV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
CPS_RV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
MRV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
RVM_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);

for c=1:size(CAL_WDW,1)
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM', ...
     'all_stn_CPS_RV', 'all_stn_EPC_RV', 'all_stn_MRV', 'all_stn_RVM')

    EPC_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
    CPS_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
    MRV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
    RVM_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
        
    for NUM_STNS = numstnstocompare
        load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox'])
        
        % Normalising
        stn_ts_nm=nan(size(stn_ts)); % nm means normalised
        stn_ts_mn=mean(stn_ts,3); 
        stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
        for n=1:NUM_TRIALS
            for m=1:NUM_STNS
               stn_ts_nm(n,m,:)= single((stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m)));
               
            end
        end
        
        % Making sure it uses positive correlated proxies (flips)
        poscorr_stn_ts = stn_ts_nm;
        mov_corr = nan(NUM_TRIALS,NUM_STNS,NUM_YRS);
        
        for n=1:NUM_TRIALS
            for m=1:NUM_STNS
                if corr(n34_ind, squeeze(stn_ts_nm(n,m,:))) < 0
                    poscorr_stn_ts(n,m,:) = -stn_ts_nm(n,m,:);
                end
                mov_corr(n,m,:) = movingCorrelation([n34_ind,squeeze(poscorr_stn_ts(n,m,:))],VAR_WDW,2);
            end
        end
        std_mov_corr = squeeze(std(mov_corr,0,2));
        EPC_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
        CPS_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
        MRV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
        RVM_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
        for n=1:NUM_TRIALS
            EPC_RV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_EPC_RV(NUM_STNS,n,VAR_WDW+1:end))));
            CPS_RV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_CPS_RV(NUM_STNS,n,VAR_WDW+1:end))));
            MRV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_MRV(NUM_STNS,n,VAR_WDW+1:end))));
            RVM_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_RVM(NUM_STNS,n,VAR_WDW+1:end))));
        end
        
    end
    toc;
    
    EPC_RV_all_grps(c,:,:) = all_stn_corr_EPC_RV;
    CPS_RV_all_grps(c,:,:) = all_stn_corr_CPS_RV;
    MRV_all_grps(c,:,:) = all_stn_corr_MRV;
    RVM_all_grps(c,:,:) = all_stn_corr_RVM;
    EPC_RV_std_mov_corr_all_grps(c,:,:) = EPC_RV_corr_std_mov_corr_ts;
    CPS_RV_std_mov_corr_all_grps(c,:,:) = CPS_RV_corr_std_mov_corr_ts;
    MRV_std_mov_corr_all_grps(c,:,:) = MRV_corr_std_mov_corr_ts;
    RVM_std_mov_corr_all_grps(c,:,:) = RVM_corr_std_mov_corr_ts; 
    
end


