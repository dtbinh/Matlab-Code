% This script will compare the RVM and MRV methods in the reconstructions

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

% 
% %% Scatterplot - Reconstruction skill vs correlation to std(corr(stn_ts,n34_ind))
% tic;
% window = 31; % The running window in years
% GROUP_NAME = 'glb_ts';
% 
% NUM_CAL_WDW = 10; clear CAL_WDW;
% overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
% for c=0:9
%     CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
% end
% 
% EPC_RV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% CPS_RV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% MRV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% RVM_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% EPC_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS);
% CPS_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS);
% MRV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS);
% RVM_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS);
% EPC_RV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% CPS_RV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% MRV_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% RVM_std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS);
% 
% for c=1%:size(CAL_WDW,1)
%     DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
%     load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
%      'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM', ...
%      'all_stn_CPS_RV', 'all_stn_EPC_RV', 'all_stn_MRV', 'all_stn_RVM')
% 
%     EPC_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%     CPS_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%     MRV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%     RVM_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%         
%     for NUM_STNS = 10
%         load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox'])
%         
%         % Normalising
%         stn_ts_nm=nan(size(stn_ts)); % nm means normalised
%         stn_ts_mn=mean(stn_ts,3); 
%         stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
%         for n=1:NUM_TRIALS
%             for m=1:NUM_STNS
%                stn_ts_nm(n,m,:)= single((stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m)));
%                
%             end
%         end
%          
%         % Making sure it uses positive correlated proxies (flips)
%         poscorr_stn_ts = stn_ts_nm;
%         mov_corr = nan(NUM_TRIALS,NUM_STNS,NUM_YRS);
%         
%         for n=1:NUM_TRIALS
%             for m=bad_recons_RVM
%                 if corr(n34_ind, squeeze(stn_ts_nm(n,m,:))) < 0
%                     poscorr_stn_ts(n,m,:) = -stn_ts_nm(n,m,:);
%                 end
%                 mov_corr(n,m,:) = movingCorrelation([n34_ind,squeeze(poscorr_stn_ts(n,m,:))],VAR_WDW,2);
%             end
%         end
%         
%         std_mov_corr = squeeze(std(mov_corr,0,2));
%         EPC_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%         CPS_RV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%         MRV_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
%         RVM_corr_std_mov_corr_ts = nan(max(numstnstocompare),NUM_TRIALS);
% %         for n=bad_recons
% %             EPC_RV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_EPC_RV(NUM_STNS,n,VAR_WDW+1:end))));
% %             CPS_RV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_CPS_RV(NUM_STNS,n,VAR_WDW+1:end))));
% %             MRV_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_MRV(NUM_STNS,n,VAR_WDW+1:end))));
% %             RVM_corr_std_mov_corr_ts(NUM_STNS,n) = corr(squeeze(std_mov_corr(n,VAR_WDW+1:end))',squeeze(double(all_stn_RVM(NUM_STNS,n,VAR_WDW+1:end))));
% %         end
%         
% 
%     end
% 
%     EPC_RV_all_grps(c,:,:,:) = all_stn_EPC_RV;
%     CPS_RV_all_grps(c,:,:,:) = all_stn_CPS_RV;
%     MRV_all_grps(c,:,:,:) = all_stn_MRV;
%     RVM_all_grps(c,:,:,:) = all_stn_RVM;
%     
%     EPC_RV_all_corr_grps(c,:,:) = all_stn_corr_EPC_RV;
%     CPS_RV_all_corr_grps(c,:,:) = all_stn_corr_CPS_RV;
%     MRV_all_corr_grps(c,:,:) = all_stn_corr_MRV;
%     RVM_all_corr_grps(c,:,:) = all_stn_corr_RVM;
%     
%     EPC_RV_std_mov_corr_all_grps(c,:,:) = EPC_RV_corr_std_mov_corr_ts;
%     CPS_RV_std_mov_corr_all_grps(c,:,:) = CPS_RV_corr_std_mov_corr_ts;
%     MRV_std_mov_corr_all_grps(c,:,:) = MRV_corr_std_mov_corr_ts;
%     RVM_std_mov_corr_all_grps(c,:,:) = RVM_corr_std_mov_corr_ts; 
%     toc;
% end
% 
% % Plotting
% subplot(2,2,1)
% scatter(EPC_RV_std_mov_corr_all_grps(:),EPC_RV_all_grps(:),1,'.')
% ylabel('Reconstruction Skill - correlation')
% xlabel('Correlation to STD time series')
% ylim([-1,1]); xlim([-1,1]);
% title('EPC\_RV');
% subplot(2,2,2)
% scatter(CPS_RV_std_mov_corr_all_grps(:),CPS_RV_all_grps(:),1,'.')
% ylabel('Reconstruction Skill - correlation')
% xlabel('Correlation to STD time series')
% ylim([-1,1]); xlim([-1,1]);
% title('CPS\_RV');
% subplot(2,2,3)
% scatter(MRV_std_mov_corr_all_grps(:),MRV_all_grps(:),1,'.')
% ylabel('Reconstruction Skill - correlation')
% xlabel('Correlation to STD time series')
% ylim([-1,1]); xlim([-1,1]);
% title('MRV');
% subplot(2,2,4)
% scatter(RVM_std_mov_corr_all_grps(:),RVM_all_grps(:),1,'.')
% ylabel('Reconstruction Skill - correlation')
% xlabel('Correlation to STD time series')
% ylim([-1,1]); xlim([-1,1]);
% title('RVM');

%% Plotting Reconstruction against std of running correlation of proxies

tic;
window = 31; % The running window in years
GROUP_NAME = 'ntrop_ts';

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

EPC_RV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,'single');
CPS_RV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,'single');
MRV_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,'single');
RVM_all_corr_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,'single');
EPC_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS,'single');
CPS_RV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS,'single');
MRV_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS,'single');
RVM_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS,'single');
std_mov_corr_all_grps = nan(size(CAL_WDW,1),max(numstnstocompare),NUM_TRIALS,NUM_YRS,'single');
bad_recons_RVM = cell(10,1); bad_recons_MRV = cell(10,1); 
bad_recons_CPS_RV = cell(10,1); bad_recons_EPC_RV = cell(10,1); 
for c=1:size(CAL_WDW,1)
    
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM', ...
     'all_stn_CPS_RV', 'all_stn_EPC_RV', 'all_stn_MRV', 'all_stn_RVM')

    EPC_RV_all_grps(c,:,:,:) = all_stn_EPC_RV;
    CPS_RV_all_grps(c,:,:,:) = all_stn_CPS_RV;
    MRV_all_grps(c,:,:,:) = all_stn_MRV;
    RVM_all_grps(c,:,:,:) = all_stn_RVM;
    EPC_RV_all_corr_grps(c,:,:) = all_stn_corr_EPC_RV;
    CPS_RV_all_corr_grps(c,:,:) = all_stn_corr_CPS_RV;
    MRV_all_corr_grps(c,:,:) = all_stn_corr_MRV;
    RVM_all_corr_grps(c,:,:) = all_stn_corr_RVM;

    for NUM_STNS = 10
        load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox'])
        
        bad_recons_RVM{c} = find(all_stn_corr_RVM(NUM_STNS,:)<prctile(squeeze(all_stn_corr_MRV(NUM_STNS,:))',5));
        bad_recons_CPS_RV{c} = find(all_stn_corr_CPS_RV(NUM_STNS,:)<prctile(squeeze(all_stn_corr_MRV(NUM_STNS,:))',5));
        bad_recons_EPC_RV{c} = find(all_stn_corr_EPC_RV(NUM_STNS,:)<prctile(squeeze(all_stn_corr_MRV(NUM_STNS,:))',5));
        bad_recons_MRV{c} = find(all_stn_corr_MRV(NUM_STNS,:)<0.5);
        % Normalising
        stn_ts_nm=nan(size(stn_ts)); poscorr_stn_ts=nan(size(stn_ts)); % nm means normalised
        stn_ts_mn=mean(stn_ts,3); 
        stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
        mov_corr = nan(NUM_TRIALS,NUM_STNS,NUM_YRS);
        for n=1:NUM_TRIALS
            for m=1:NUM_STNS
               stn_ts_nm(n,m,:)= single((stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m)));
               if corr(n34_ind, squeeze(stn_ts_nm(n,m,:))) < 0
                    poscorr_stn_ts(n,m,:) = -stn_ts_nm(n,m,:);
               else
                   poscorr_stn_ts(n,m,:) = stn_ts_nm(n,m,:);
               end
               mov_corr(n,m,:) = movingCorrelation([n34_ind,squeeze(poscorr_stn_ts(n,m,:))],VAR_WDW,2);
            end
        end
        
        std_mov_corr_all_grps(c,NUM_STNS,:,:) = squeeze(std(mov_corr,0,2));
        toc;
    end
end

save('DataFiles/rvm_vs_mrv_ntrop.mat','std_mov_corr_all_grps','bad_recons_CPS_RV','bad_recons_EPC_RV','bad_recons_RVM','bad_recons_MRV',...
       'EPC_RV_all_grps','CPS_RV_all_grps','MRV_all_grps','RVM_all_grps','EPC_RV_all_corr_grps','CPS_RV_all_corr_grps',...
       'MRV_all_corr_grps','RVM_all_corr_grps','window','GROUP_NAME');
%% Plotting
c=10;
i=5;

s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_MRV{c}(i),(31:end)));
rec_skill = squeeze(MRV_all_corr_grps(c,NUM_STNS,bad_recons_MRV{c}(i)));
r = squeeze(MRV_all_grps(c,NUM_STNS,bad_recons_MRV{c}(i),(16:end-15)));
% scatter(std_mov_corr(31:end),recon(16:end-15))
subplot(4,1,1)
scatter(s(:),r(:),'.');
% values = hist3([s(:) r(:)],{0:0.01:1 , 0:0.01:3})';
% pcolor(0:0.01:1,0:0.01:3,values); shading flat; colormap(flipud(gray)); colorbar;
xlabel('std of running correlation of proxies to ENSO')
ylabel('Reconstructed Nino3.4 SST RV')
xlim([0 0.5]); ylim([0 2.5]); grid on;
[b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
hold on; plot(linspace(0,0.5,100),b(2)+linspace(0,0.5,100)*b(1),'k'); hold off
text(0.4,1.25,['R^2 = ',num2str(thestats(1))])
text(0.4,1,['P-value = ',num2str(thestats(4))])
text(0.4,1.75,['y=',num2str(b(2)),'+',num2str(b(1)),' x'])
text(0.4,1.5,['r(MRV,n34) = ',num2str(RVM_all_corr_grps(c,NUM_STNS,bad_recons_RVM{c}(i)))])
title('MRV')
% std_mov_corr = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_RVM(i),:));
% recon = squeeze(RVM_all_grps(NUM_STNS,bad_recons_RVM(i),:));
% plotyy(1:469,recon(16:end-15),1:469,std_mov_corr(31:end));
% corr(std_mov_corr(31:end),recon(16:end-15))
s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_RVM{c}(i),(31:end)));
rec_skill = squeeze(RVM_all_corr_grps(c,NUM_STNS,bad_recons_RVM{c}(i)));
r = squeeze(RVM_all_grps(c,NUM_STNS,bad_recons_RVM{c}(i),(16:end-15)));
% scatter(std_mov_corr(31:end),recon(16:end-15))
subplot(4,1,2)
scatter(s(:),r(:),'.');
% values = hist3([s(:) r(:)],{0:0.01:1 , 0:0.01:3})';
% pcolor(0:0.01:1,0:0.01:3,values); shading flat; colormap(flipud(gray)); colorbar;
xlabel('std of running correlation of proxies to ENSO')
ylabel('Reconstructed Nino3.4 SST RV')
xlim([0 0.5]); ylim([0 2.5]); grid on;
[b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
hold on; plot(linspace(0,0.5,100),b(2)+linspace(0,0.5,100)*b(1),'k'); hold off
text(0.4,1.25,['R^2 = ',num2str(thestats(1))])
text(0.4,1,['P-value = ',num2str(thestats(4))])
text(0.4,1.75,['y=',num2str(b(2)),'+',num2str(b(1)),' x'])
text(0.4,1.5,['r(RVM,n34) = ',num2str(RVM_all_corr_grps(c,NUM_STNS,bad_recons_RVM{c}(i)))])
title('RVM')

s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i),(31:end)));
rec_skill = squeeze(CPS_RV_all_corr_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i)));
r = squeeze(CPS_RV_all_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i),(16:end-15)));
% scatter(std_mov_corr(31:end),recon(16:end-15))
subplot(4,1,3)
scatter(s(:),r(:),'.');
% values = hist3([s(:) r(:)],{0:0.01:1 , 0:0.01:3})';
% pcolor(0:0.01:1,0:0.01:3,values); shading flat; colormap(flipud(gray)); colorbar;
xlabel('std of running correlation of proxies to ENSO')
ylabel('Reconstructed Nino3.4 SST RV')
xlim([0 0.5]); ylim([0 2.5]); grid on;
[b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
hold on; plot(linspace(0,0.5,100),b(2)+linspace(0,0.5,100)*b(1),'k'); hold off
text(0.4,1.25,['R^2 = ',num2str(thestats(1))])
text(0.4,1.75,['y=',num2str(b(2)),'+',num2str(b(1)),' x'])
text(0.4,1,['P-value = ',num2str(thestats(4))])
text(0.4,1.5,['r(CPS\_RV,n34) = ',num2str(RVM_all_corr_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i)))])
title('CPS\_RV')

s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i),(31:end)));
rec_skill = squeeze(EPC_RV_all_corr_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i)));
r = squeeze(EPC_RV_all_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i),(16:end-15)));
% scatter(std_mov_corr(31:end),recon(16:end-15))
subplot(4,1,4)
scatter(s(:),r(:),'.');
% values = hist3([s(:) r(:)],{0:0.01:1 , 0:0.01:3})';
% pcolor(0:0.01:1,0:0.01:3,values); shading flat; colormap(flipud(gray)); colorbar;
xlabel('std of running correlation of proxies to ENSO')
ylabel('Reconstructed Nino3.4 SST RV')
xlim([0 0.5]); ylim([0 2.5]); grid on;
[b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
hold on; plot(linspace(0,0.5,100),b(2)+linspace(0,0.5,100)*b(1),'k'); hold off
text(0.4,1.25,['R^2 = ',num2str(thestats(1))])
text(0.4,1.75,['y=',num2str(b(2)),'+',num2str(b(1)),' x'])
text(0.4,1,['P-value = ',num2str(thestats(4))])
text(0.4,1.5,['r(EPC\_RV,n34) = ',num2str(RVM_all_corr_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i)))])
title('EPC\_RV')

%% Getting alot of R^2 and regression values


as_RVM=[]; ar_RVM=[];  ab_RVM=[];  ap_RVM=[]; ac_RVM=[];
as_MRV=[]; ar_MRV=[];  ab_MRV=[];  ap_MRV=[]; ac_MRV=[];
as_CPS_RV=[]; ar_CPS_RV=[];  ab_CPS_RV=[];  ap_CPS_RV=[]; ac_CPS_RV=[];
as_EPC_RV=[]; ar_EPC_RV=[];  ab_EPC_RV=[];  ap_EPC_RV=[]; ac_EPC_RV=[];

for c=1:10

    all_skill_RVM=nan(length(bad_recons_RVM{c}(:)),1,'single');
    allrsq_RVM=nan(length(bad_recons_RVM{c}(:)),1,'single');
    allb_RVM=nan(length(bad_recons_RVM{c}(:)),2,'single');
    allp_RVM=nan(length(bad_recons_RVM{c}(:)),1,'single');
    allc_RVM=nan(length(bad_recons_RVM{c}(:)),1,'single'); % Correlation to std
    
    all_skill_MRV=nan(length(bad_recons_MRV{c}(:)),1,'single');
    allrsq_MRV=nan(length(bad_recons_MRV{c}(:)),1,'single');
    allb_MRV=nan(length(bad_recons_MRV{c}(:)),2,'single');
    allp_MRV=nan(length(bad_recons_MRV{c}(:)),1,'single');
    allc_MRV=nan(length(bad_recons_MRV{c}(:)),1,'single');
    
    all_skill_CPS_RV=nan(length(bad_recons_CPS_RV{c}(:)),1,'single');
    allrsq_CPS_RV=nan(length(bad_recons_CPS_RV{c}(:)),1,'single');
    allb_CPS_RV=nan(length(bad_recons_CPS_RV{c}(:)),2,'single');
    allp_CPS_RV=nan(length(bad_recons_CPS_RV{c}(:)),1,'single');
    allc_CPS_RV=nan(length(bad_recons_CPS_RV{c}(:)),1,'single');
    
    all_skill_EPC_RV=nan(length(bad_recons_EPC_RV{c}(:)),1,'single');
    allrsq_EPC_RV=nan(length(bad_recons_EPC_RV{c}(:)),1,'single');
    allb_EPC_RV=nan(length(bad_recons_EPC_RV{c}(:)),2,'single');
    allp_EPC_RV=nan(length(bad_recons_EPC_RV{c}(:)),1,'single');
    allc_EPC_RV=nan(length(bad_recons_EPC_RV{c}(:)),1,'single');

    for i=1:length(bad_recons_MRV{c}(:));
    s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_MRV{c}(i),(31:end)));
    r = squeeze(MRV_all_grps(c,NUM_STNS,bad_recons_MRV{c}(i),(16:end-15)));
    all_skill_MRV(i) = squeeze(MRV_all_corr_grps(c,NUM_STNS,bad_recons_MRV{c}(i)));
    [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
    allrsq_MRV(i) = thestats(1);
    allp_MRV(i) = thestats(4);
    allb_MRV(i,:) = b;
    allc_MRV(i) = corr(r(:),s(:));
    end

    for i=1:length(bad_recons_RVM{c}(:));
    s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_RVM{c}(i),(31:end)));
    r = squeeze(RVM_all_grps(c,NUM_STNS,bad_recons_RVM{c}(i),(16:end-15)));
    all_skill_RVM(i) = squeeze(RVM_all_corr_grps(c,NUM_STNS,bad_recons_RVM{c}(i)));
    [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
    allrsq_RVM(i) = thestats(1);
    allp_RVM(i) = thestats(4);
    allb_RVM(i,:) = b;
    allc_RVM(i) = corr(r(:),s(:));
    end

    for i=1:length(bad_recons_CPS_RV{c}(:));
    s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i),(31:end)));
    r = squeeze(CPS_RV_all_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i),(16:end-15)));
    all_skill_CPS_RV(i) = squeeze(CPS_RV_all_corr_grps(c,NUM_STNS,bad_recons_CPS_RV{c}(i)));
    [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
    allrsq_CPS_RV(i) = thestats(1);
    allb_CPS_RV(i,:) = b;
    allp_CPS_RV(i) = thestats(4);
    allc_CPS_RV(i) = corr(r(:),s(:));
    end

    for i=1:length(bad_recons_EPC_RV{c}(:));
    s = squeeze(std_mov_corr_all_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i),(31:end)));
    r = squeeze(EPC_RV_all_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i),(16:end-15)));
    all_skill_EPC_RV(i) = squeeze(EPC_RV_all_corr_grps(c,NUM_STNS,bad_recons_EPC_RV{c}(i)));
    [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
    allrsq_EPC_RV(i) = thestats(1);
    allb_EPC_RV(i,:) = b;
    allp_EPC_RV(i) = thestats(4);
    allc_EPC_RV(i) = corr(r(:),s(:));
    end


    as_MRV = cat(1,as_MRV,all_skill_MRV);
    ar_MRV = cat(1,ar_MRV,allrsq_MRV);
    ab_MRV = cat(1,ab_MRV,allb_MRV);
    ap_MRV = cat(1,ap_MRV,allp_MRV);
    ac_MRV = cat(1,ac_MRV,allc_MRV);
    
    as_RVM = cat(1,as_RVM,all_skill_RVM);
    ar_RVM = cat(1,ar_RVM,allrsq_RVM);
    ab_RVM = cat(1,ab_RVM,allb_RVM);
    ap_RVM = cat(1,ap_RVM,allp_RVM);
    ac_RVM = cat(1,ac_RVM,allc_RVM);
    
    as_CPS_RV = cat(1,as_CPS_RV,all_skill_CPS_RV);
    ar_CPS_RV = cat(1,ar_CPS_RV,allrsq_CPS_RV);
    ab_CPS_RV = cat(1,ab_CPS_RV,allb_CPS_RV);
    ap_CPS_RV = cat(1,ap_CPS_RV,allp_CPS_RV);
    ac_CPS_RV = cat(1,ac_CPS_RV,allc_CPS_RV);
    
    as_EPC_RV = cat(1,as_EPC_RV,all_skill_EPC_RV);
    ar_EPC_RV = cat(1,ar_EPC_RV,allrsq_EPC_RV);
    ab_EPC_RV = cat(1,ab_EPC_RV,allb_EPC_RV);
    ap_EPC_RV = cat(1,ap_EPC_RV,allp_EPC_RV);
    ac_EPC_RV = cat(1,ac_EPC_RV,allc_EPC_RV);

end

% Finding when reg is significant

SIGREG_CUTOFF = -1; % For regression coefficient
SIGRSQ_CUTOFF = 0.2;
SIGP_CUTOFF = 0.05;

sig_MRV = find(ap_MRV<SIGP_CUTOFF);
sig_RVM = find(ap_RVM<SIGP_CUTOFF);
sig_CPS_RV = find(ap_CPS_RV<SIGP_CUTOFF);
sig_EPC_RV = find(ap_EPC_RV<SIGP_CUTOFF);

% Corr to std vs Reg coef
clf
subplot(2,2,1)
plot(squeeze(ab_MRV(:,1)),ac_MRV(:),'b.'); hold on;
plot(squeeze(ab_MRV(sig_MRV,1)),ac_MRV(sig_MRV),'g.');
xlabel('regression coef')
ylabel('Corr value')
title('MRV','FontSize',14)
xlim([-7 5]); ylim([-0.8 0.8]); grid on
legend('All poor reconstructions','Significant reconstructions');
h=text(-1,0.4,[num2str(length(sig_MRV)/length(ap_MRV)*100,'%3.1f'),'% sig']);
set(h,'FontSize',14)

subplot(2,2,2)
plot(squeeze(ab_RVM(:,1)),ac_RVM(:),'b.'); hold on;
plot(squeeze(ab_RVM(sig_RVM,1)),ac_RVM(sig_RVM),'g.')
xlabel('regression coef')
ylabel('Corr value')
title('RVM','FontSize',14)
xlim([-7 5]); ylim([-0.8 0.8]); grid on
h=text(-1,0.4,[num2str(length(sig_RVM)/length(ap_RVM)*100,'%3.1f'),'% sig']);
set(h,'FontSize',14)

subplot(2,2,3)
plot(squeeze(ab_CPS_RV(:,1)),ac_CPS_RV(:),'b.'); hold on;
plot(squeeze(ab_CPS_RV(sig_CPS_RV,1)),ac_CPS_RV(sig_CPS_RV),'g.')
xlabel('regression coef')
ylabel('Corr value')
title('CPS\_RV','FontSize',14)
xlim([-7 5]); ylim([-0.8 0.8]); grid on
h=text(-1,0.4,[num2str(length(sig_CPS_RV)/length(ap_CPS_RV)*100,'%3.1f'),'% sig']);
set(h,'FontSize',14)

subplot(2,2,4)
plot(squeeze(ab_EPC_RV(:,1)),ac_EPC_RV(:),'b.'); hold on;
plot(squeeze(ab_EPC_RV(sig_EPC_RV,1)),ac_EPC_RV(sig_EPC_RV),'g.')
xlabel('regression coef')
ylabel('Corr value')
title('EPC\_RV','FontSize',14)
xlim([-7 5]); ylim([-0.8 0.8]); grid on
h=text(-1,0.4,[num2str(length(sig_EPC_RV)/length(ap_EPC_RV)*100,'%3.1f'),'% sig']);
set(h,'FontSize',14)

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/corr(rec,std(prox))_vs_regcoef_badrecons.jpg'])

% Reg coef vs recon skill
clf;
subplot(4,1,1)
plot(squeeze(ab_MRV(sig_MRV,1)),as_MRV(sig_MRV),'g.')
lsline
xlabel('regression coef')
ylabel('skill')
title('MRV')
xlim([-7 3]); ylim([-0.1 1]);

subplot(4,1,2)
plot(squeeze(ab_RVM(sig_RVM,1)),as_RVM(sig_RVM),'g.')
lsline
xlabel('regression coef')
ylabel('skill')
title('RVM')
xlim([-7 3]); ylim([-0.1 1]);

subplot(4,1,3)
plot(squeeze(ab_CPS_RV(sig_CPS_RV,1)),as_CPS_RV(sig_CPS_RV),'g.')
lsline
xlabel('regression coef')
ylabel('skill')
title('CPS\_RV')
xlim([-7 3]); ylim([-0.1 1]);

subplot(4,1,4)
plot(squeeze(ab_EPC_RV(sig_EPC_RV,1)),as_EPC_RV(sig_EPC_RV),'g.')
lsline
xlabel('regression coef')
ylabel('skill')
title('EPC\_RV')
xlim([-7 3]); ylim([-0.1 1]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/regcoef_vs_skill_badrecons.jpg'])

% recon skill vs Corr to std
clf;
subplot(4,1,1)
plot(ac_MRV(sig_MRV),as_MRV(sig_MRV),'g.')
lsline
xlabel('Corr')
ylabel('skill')
title('MRV')
xlim([-0.8 0.8]); ylim([-0.1 0.6]);

subplot(4,1,2)
plot(ac_RVM(sig_RVM),as_RVM(sig_RVM),'g.')
lsline
xlabel('Corr')
ylabel('skill')
title('RVM')
xlim([-0.8 0.8]); ylim([-0.1 0.6]);

subplot(4,1,3)
plot(ac_CPS_RV(sig_CPS_RV),as_CPS_RV(sig_CPS_RV),'g.')
lsline
xlabel('Corr')
ylabel('skill')
title('CPS\_RV')
xlim([-0.8 0.8]); ylim([-0.1 0.6]);

subplot(4,1,4)
plot(ac_EPC_RV(sig_EPC_RV),as_EPC_RV(sig_EPC_RV),'g.')
lsline
xlabel('Corr')
ylabel('skill')
title('EPC\_RV')
xlim([-0.8 0.8]); ylim([-0.1 0.6]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/corr(rec,std(prox))_vs_skill_badrecons.jpg'])