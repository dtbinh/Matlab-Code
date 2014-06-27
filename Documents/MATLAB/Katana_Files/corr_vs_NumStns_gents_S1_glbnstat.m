% This script will produce a plot of the percentiles of correlation
% corresponding to the number of stations. This is mainly a combination of
% code in recons.m and station_select.m, but modified for many station
% output. This is Series 1 of experiments.
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
% window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];

%% Beggining of Loop

DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/glb_ts_nstat_sigpcd'] ; mkdir(DIR_NAME);
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']);
load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat']);
% load DataFiles/runcorr_eofs.mat

NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
% lsfrac=nc_varget('DataFiles/sftlf_A1.static.nc','sftlf'); lsfrac(isnan(lsfrac)) = 0 ;
% Random selection boundaries

% S_lat = -50; N_lat = 15; W_lon = 280; E_lon = 330; % South America
% S_lat = -40; N_lat = -12; W_lon = 110; E_lon = 160; % Australia
% S_lat = 10; N_lat = 50; W_lon = 235; E_lon = 280; % USA
S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360; % Global
% S_lat = -60; N_lat = 60; W_lon = 0; E_lon = 360; % Nonpolar
% S_lat = -10; N_lat = 10; W_lon = 100; E_lon = 300; % Tropical
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

% Selection from areas with absolute correlation over a certain threshold
MIN_COR = 0.3;
% CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];
CAL_WDW = 1:499;

% Changing correlation in correlation matrix so that certain areas are not
% picked up in the station selection

% [~,S_bd]= min(abs(lat--10));
% [~,N_bd]= min(abs(lat-10));
% [~,W_bd]= min(abs(lon-100));
% [~,E_bd]= min(abs(lon-300));
% 
% corr_ts(S_bd:N_bd,W_bd:E_bd) = 0;
STN_MAX = 70;
indice_pool_num = zeros(STN_MAX,1);
%% Selecting stations

corr_ts = nan(size(ats,2),size(ats,3),'single');
for i=S_bound:N_bound
    for j=W_bound:E_bound
        corr_ts(i,j) = corr(n34_ind(CAL_WDW),ats(CAL_WDW,i,j));
    end
end

mkdir([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(1)),'-',num2str(CAL_WDW(end))])

for NUM_STNS = 3:STN_MAX
    
    stn_ts = nan(NUM_TRIALS,NUM_STNS, NUM_YRS,'single');
    stn_pr = nan(NUM_TRIALS,NUM_STNS, NUM_YRS,'single');
    stn_lat = nan(NUM_TRIALS,NUM_STNS);
    stn_lon = nan(NUM_TRIALS,NUM_STNS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conditions - Replace 1 with desired conditions

    indice_pool = find(abs(corr_ts)>MIN_COR & ...
                       nonstat_tsmap > ceil(0.1*(NUM_YRS-window)) & ...
                       squeeze(abs(mean(ts_pc(1,(window+2):end,:,:)-ts_pc(2,(window+2):end,:,:),2))) > 0.3 & ...
                       1 & ...
                       1                                                     );

    %  - squeeze(abs(mean(ts_pc(1,33:end,:,:)-ts_pc(2,33:end,:,:),2))) > 0.3
    %  - nonstat_tsmap > 50
    %  - lsfrac > 0
    %  - squeeze(eof_ts_fm(1,:,:)) > 0.01
    %  - squeeze(eof_ts_fm(2,:,:)) < -0.01
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:NUM_TRIALS

        [stn_lat(m,:),stn_lon(m,:)] = ind2sub(size(corr_ts),indice_pool(randperm(length(indice_pool),NUM_STNS)));

        for n=1:NUM_STNS

            stn_ts(m,n,:) = single(ats(:,stn_lat(m,n),stn_lon(m,n)));
            stn_pr(m,n,:) = single(apr(:,stn_lat(m,n),stn_lon(m,n)));

        end
    end
    
    save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(1)),'-',num2str(CAL_WDW(end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
             'stn_lat','stn_lon','stn_ts','indice_pool','corr_ts','window');

end

% Writing README file

fid = fopen([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(1)),'-',num2str(CAL_WDW(end)),'/README.txt'], 'w+');

fprintf(fid,'%% This file contains the parameters that were used with corr_vs_NumStns_gen.m ');
fprintf(fid,['when generating the files in this folder(',DIR_NAME,').\n\n']);
fprintf(fid,'Boundaries of the box where the stations are contained.\n');
fprintf(fid,['S_lat = ',num2str(S_lat),'\n']);
fprintf(fid,['N_lat = ',num2str(N_lat),'\n']);
fprintf(fid,['W_lon = ',num2str(W_lon),'\n']);
fprintf(fid,['E_lon = ',num2str(E_lon),'\n\n']);
fprintf(fid,'Selection from areas with absolute correlation over a certain threshold.\n');
fprintf(fid,['MIN_COR = ',num2str(MIN_COR),'\n']);
fprintf(fid,['CAL_WDW = ',num2str(CAL_WDW(1)),':',num2str(CAL_WDW(end)),'\n\n']);
fprintf(fid,'Data is also using temperature only\n');
% fprintf(fid,'Data is also using precipitation only\n');
fprintf(fid,'This file was produced using the UNSW Katana Computational Cluster.\n\n');
fprintf(fid,'Station selection conditions:\n');
fprintf(fid,'abs(corr_ts)>MIN_COR\n');
fprintf(fid,'squeeze(abs(mean(ts_pc(1,(window+2):end,:,:)-ts_pc(2,(window+2):end,:,:),2))) > 0.3\n');
fprintf(fid,'nonstat_tsmap > ceil(0.1*(NUM_YRS-window))\n');
% fprintf(fid,'lsfrac > 0\n');
% fprintf(fid,'Tropical regions have not been included\n');
% fprintf(fid,'eof_ts_fm2 > 0.01\n');
fclose(fid);

% %% Data Processing and Reconstruction Methods
% tic;
% STN_NUM_RG = [3:STN_MAX];
% all_stn_corr_MRV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_rmse_MRV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_rmse_EPC = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_corr_EPC = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_rmse_EPC_RV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_corr_EPC_RV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_rmse_CPS_RV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_corr_CPS_RV = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_rmse_CPS = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% all_stn_corr_CPS = nan(max(STN_NUM_RG), NUM_TRIALS,'single');
% 
% for NUM_STNS = STN_NUM_RG
%     
%     % Loading and normalising proxies
%     load([DIR_NAME,'/',num2str(NUM_STNS),'stns_1000prox.mat']);
%     stn_ts_mn=mean(stn_ts,3); 
%     stn_ts_std=squeeze(std(permute(stn_ts,[3,1,2])));
%     for n=1:NUM_TRIALS
%         for m=1:NUM_STNS
%            stn_ts(n,m,:)= (stn_ts(n,m,:)-stn_ts_mn(n,m))./(stn_ts_std(n,m));
%         end
%     end
%     
%     
%     % McGregor et al 2013 method of Median Running Variances
%     stn_MRV=nan(NUM_TRIALS, NUM_YRS);
%     for n=1:NUM_TRIALS
%         stn_movvar=nan(NUM_STNS,NUM_YRS);
%         for m=1:NUM_STNS
%             stn_movvar(m,:) = movingvar(squeeze(stn_ts(n,m,:)),VAR_WDW);
%         end
%         stn_MRV(n,:) = median(stn_movvar);
%     end
%     % Skill Evaluation
%     stn_corr_MRV = nan(NUM_TRIALS,1);
%     stn_rmse_MRV = nan(NUM_TRIALS,1);
%     for n=1:NUM_TRIALS
%         stn_corr_MRV(n) = corr(squeeze(stn_MRV(n,RV_WDW))',n34_ind_RV(RV_WDW));
%         stn_rmse_MRV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_MRV(n,RV_WDW)').^2));
%     end
%     all_stn_corr_MRV(NUM_STNS,:) = stn_corr_MRV;
%     all_stn_rmse_MRV(NUM_STNS,:) = stn_rmse_MRV;
%     
%     
%     % Braganza et al 2009 method of EOF construction
%     NUM_OF_EOFS = 3;
%     eof_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_STNS);
%     PC_stn = nan(NUM_TRIALS,NUM_OF_EOFS,NUM_YRS);
%     expvar_stn = nan(NUM_TRIALS,NUM_OF_EOFS);
%     for n=1:NUM_TRIALS
%         [eof_stn(n,:,:),PC_stn(n,:,:),expvar_stn(n,:)] = caleof(squeeze(stn_ts(n,:,:))', NUM_OF_EOFS, 1);
%     end
%     % Flipping EOFs where necessary so PC is similar to PC_stn(1,1,:)
%     for n=2:NUM_TRIALS
%         if corr(squeeze(PC_stn(1,1,:)), squeeze(PC_stn(n,1,:))) < 0
%             PC_stn(n,1,:) = -PC_stn(n,1,:);
%             eof_stn(n,1,:) = -eof_stn(n,1,:);
%         end
%     end
%     % Normalising (it already has mean 0)
%     stn_EPC = nan(NUM_TRIALS, NUM_YRS);
%     stn_EPC_RV = nan(NUM_TRIALS, NUM_YRS);
%     for n=1:NUM_TRIALS
%         stn_EPC(n,:) = squeeze(PC_stn(n,1,:))./std(squeeze(PC_stn(n,1,:))');
%         stn_EPC_RV(n,:) = movingvar(stn_EPC(n,:)',VAR_WDW);
%     end
%     % Skill Evaluation
%     stn_corr_EPC = nan(NUM_TRIALS,1);
%     stn_rmse_EPC = nan(NUM_TRIALS,1);
%     stn_corr_EPC_RV = nan(NUM_TRIALS,1);
%     stn_rmse_EPC_RV = nan(NUM_TRIALS,1);
%     for n=1:NUM_TRIALS
%         stn_corr_EPC(n) = abs(corr(squeeze(stn_EPC(n,:))',n34_ind));
%         stn_rmse_EPC(n) = sqrt(mean((n34_ind-stn_EPC(n,:)').^2));
%         stn_corr_EPC_RV(n) = abs(corr(squeeze(stn_EPC_RV(n,RV_WDW))',n34_ind_RV(RV_WDW)));
%         stn_rmse_EPC_RV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_EPC_RV(n,RV_WDW)').^2));
%     end
%     all_stn_corr_EPC_RV(NUM_STNS,:) = stn_corr_EPC_RV;
%     all_stn_rmse_EPC_RV(NUM_STNS,:) = stn_rmse_EPC_RV;
%     all_stn_corr_EPC(NUM_STNS,:) = stn_corr_EPC;
%     all_stn_rmse_EPC(NUM_STNS,:) = stn_rmse_EPC;
%     
%     
%     % Esper et al 2005 CPS Method 
%     CAL_WDW = [1:50];
%     stn_CPS = nan(NUM_TRIALS, NUM_YRS);
%     stn_CPS_RV = nan(NUM_TRIALS, NUM_YRS);
%     for n=1:NUM_TRIALS
%         corr_matrix = corr(n34_ind(CAL_WDW)*ones(1,10), squeeze(stn_ts(n,:,CAL_WDW))');
%         stn_CPS(n,:) = corr_matrix(1,:)*squeeze(stn_ts(n,:,:));
%     end
%     % Normalising (it already has mean ~0)
%     for n=1:NUM_TRIALS
%         stn_CPS(n,:) = squeeze(stn_CPS(n,:))./std(squeeze(stn_CPS(n,:))');
%         stn_CPS_RV(n,:) = movingvar(squeeze(stn_CPS(n,:))',VAR_WDW);
%     end
%     % Skill Evaluation
%     stn_corr_CPS = nan(NUM_TRIALS,1);
%     stn_rmse_CPS = nan(NUM_TRIALS,1);
%     stn_corr_CPS_RV = nan(NUM_TRIALS,1);
%     stn_rmse_CPS_RV = nan(NUM_TRIALS,1);
%     for n=1:NUM_TRIALS
%         stn_corr_CPS(n) = corr(squeeze(stn_CPS(n,:))',n34_ind);
%         stn_rmse_CPS(n) = sqrt(mean((n34_ind-stn_CPS(n,:)').^2));
%         stn_corr_CPS_RV(n) = corr(squeeze(stn_CPS_RV(n,RV_WDW))',n34_ind_RV(RV_WDW));
%         stn_rmse_CPS_RV(n) = sqrt(mean((n34_ind_RV(RV_WDW)-stn_CPS_RV(n,RV_WDW)').^2));
%     end
%     all_stn_corr_CPS_RV(NUM_STNS,:) = stn_corr_CPS_RV;
%     all_stn_rmse_CPS_RV(NUM_STNS,:) = stn_rmse_CPS_RV;
%     all_stn_corr_CPS(NUM_STNS,:) = stn_corr_CPS;
%     all_stn_rmse_CPS(NUM_STNS,:) = stn_rmse_CPS;
%     toc;
% end % Took 40 min for 70 
% 
% save([DIR_NAME,'/method_stats.mat'],'all_stn_corr_MRV','all_stn_rmse_MRV',...
%      'all_stn_corr_EPC', 'all_stn_rmse_EPC', ...
%      'all_stn_corr_EPC_RV', 'all_stn_rmse_EPC_RV', ...
%      'all_stn_corr_CPS', 'all_stn_rmse_CPS', ...
%      'all_stn_corr_CPS_RV', 'all_stn_rmse_CPS_RV' , 'indice_pool_num','window');

% %% Making the plot
% 
% % Percentiles
% load([DIR_NAME,'/method_stats.mat'])
% point_NUM_STNS = [1:STN_MAX]'*ones(1,1000);
% % scatter(point_NUM_STNS(:),all_stn_corr_MRV(:),'.')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 28 19]); %x_width=28cm y_width=19cm
% % Multi Color composition
% subplot(1,3,1)
% quan_plot = quantile(all_stn_corr_MRV',[0.01:0.01:1]);
% qs = [1 5; 6 25; 26 50; 51 75; 76 95; 96 100];
% for n=1:6
%     scatter(reshape(point_NUM_STNS(:,qs(n,1):qs(n,2)),[(qs(n,2)-qs(n,1)+1)*STN_MAX 1]), ...
%             reshape(quan_plot(qs(n,1):qs(n,2),:)',[(qs(n,2)-qs(n,1)+1)*STN_MAX,1]),'.');
% hold on
% end
% xlim([0 STN_MAX]);
% ylim([0 1]);
% title('all\_stn\_corr\_MRV');
% xlabel('Number of Stations in Reconstruction')
% ylabel('Correlation')
% hold off; grid on;
% 
% subplot(1,3,2)
% quan_plot = quantile(all_stn_corr_EPC_RV',[0.01:0.01:1]);
% qs = [1 5; 6 25; 26 50; 51 75; 76 95; 96 100];
% for n=1:6
%     scatter(reshape(point_NUM_STNS(:,qs(n,1):qs(n,2)),[(qs(n,2)-qs(n,1)+1)*STN_MAX 1]), ...
%             reshape(quan_plot(qs(n,1):qs(n,2),:)',[(qs(n,2)-qs(n,1)+1)*STN_MAX,1]),'.');
% hold on
% end
% title('all\_stn\_corr\_EPC\_RV');
% xlim([0 STN_MAX]);
% ylim([0 1]);
% xlabel('Number of Stations in Reconstruction')
% hold off; grid on;
% 
% subplot(1,3,3)
% quan_plot = quantile(all_stn_corr_CPS_RV',[0.01:0.01:1]);
% qs = [1 5; 6 25; 26 50; 51 75; 76 95; 96 100];
% for n=1:6
%     scatter(reshape(point_NUM_STNS(:,qs(n,1):qs(n,2)),[(qs(n,2)-qs(n,1)+1)*STN_MAX 1]), ...
%             reshape(quan_plot(qs(n,1):qs(n,2),:)',[(qs(n,2)-qs(n,1)+1)*STN_MAX,1]),'.');
% hold on
% end
% title('all\_stn\_corr\_CPS\_RV');
% xlim([0 STN_MAX]);
% xlabel('Number of Stations in Reconstruction')
% hold off; grid on;
% suptitle(strrep(DIR_NAME,'_','\_'));
% legend('1-5th','6-25th','26-50th','51-75th','76-95th','96-100th','Location','SouthEast')
% ylim([0 1])
% saveas(gcf,['Plots/corr_vs_NUM_STNS_',DIR_NAME(15:end),'.jpg'])
% 
% % This takes 40 min.
% 
% 
% %% Plotting locations of stations
% ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
% lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
% load Pseudoproxies/glb/CalWdw:101-150/10stns_1000prox.mat
% plotworld; hold on; scatter(lon(stn_lon(:)),lat(stn_lat(:)),'.');
% % load Pseudoproxies/glb_neof2/10stns_1000prox.mat
% % scatter(lon(stn_lon(:)),lat(stn_lat(:)),'.');
% % load Pseudoproxies/glb_land/3stns_1000prox.mat
% % scatter(lon(stn_lon(:)),lat(stn_lat(:)),'.');
% xlim([0 360]); ylim([-90 90]);
