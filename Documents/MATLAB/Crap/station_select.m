% This script will select specific areas for the use for proxy
% reconstructions

%% Setup
ts_file = 'ts_A1.nc'; pr_file = 'pr_A1.nc';

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

window = 31; % The running window in years

%% Selecting stations


load runcorr.mat
load nonstat_map.mat

NUM_STNS = 10; NUM_SYNRUNS = 1000; NUM_YRS = 499; NUM_TRIALS = 1000;
DIR_NAME = 'Pseudoproxies_highlat' ; mkdir(DIR_NAME);
lsfrac=nc_varget('sftlf_A1.static.nc','sftlf'); lsfrac(isnan(lsfrac)) = 0 ;
% stn_synts = nan(NUM_STNS,NUM_SYNRUNS,NUM_YRS);
% stn_synpr = nan(NUM_STNS,NUM_SYNRUNS,NUM_YRS);
stn_ts = nan(NUM_STNS, NUM_YRS);
stn_pr = nan(NUM_STNS, NUM_YRS);
% stn_ts_pc = nan(NUM_STNS, 2, NUM_YRS);
% stn_pr_pc = nan(NUM_STNS, 2, NUM_YRS);

% Random selection boundaries

S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360; % Global
% S_lat = -60; N_lat = 60; W_lon = 0; E_lon = 360; % Nonpolar
% S_lat = -10; N_lat = 10; W_lon = 100; E_lon = 300; % Tropical
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

% Selection from areas with absolute correlation over a certain threshold
MIN_COR = 0.3;
COR_WDW = [1:499]; % The specific years

corr_ts = nan(size(ats,2),size(ats,3));
for i=S_bound:N_bound
    for j=W_bound:E_bound
        corr_ts(i,j) = corr(n34_ind(COR_WDW),ats(COR_WDW,i,j));
    end
end

% Changing correlation in correlation matrix so that certain areas are not
% picked up in the station selection

[~,S_bd]= min(abs(lat--10));
[~,N_bd]= min(abs(lat-10));
[~,W_bd]= min(abs(lon-100));
[~,E_bd]= min(abs(lon-300));

corr_ts(S_bd:N_bd,W_bd:E_bd) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conditions - Replace 1 with desired conditions

indice_pool = find(abs(corr_ts)>MIN_COR & ...
                   1 & ...
                   1 & ...
                   1 & ...
                   1                                                     );

%  - squeeze(abs(mean(ts_pc(1,33:end,:,:)-ts_pc(2,33:end,:,:),2))) > 0.3
%  - nonstat_tsmap > 50
%  - lsfrac > 0.75

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:NUM_TRIALS

    [stn_lat,stn_lon] = ind2sub(size(corr_ts),indice_pool(randperm(length(indice_pool),NUM_STNS)));
    
    for n=1:NUM_STNS

        stn_ts(n,:) = ats(:,stn_lat(n),stn_lon(n));
        stn_pr(n,:) = apr(:,stn_lat(n),stn_lon(n));
%         load(['Synth_corr_',num2str(window),'yrwin/',num2str(lon(stn_lon(n))),'E',num2str(lat(stn_lat(n))),'N_syncorr.mat'])
%         stn_synts(n,:,:) = spot_ts;
%         stn_synpr(n,:,:) = spot_pr;
%         stn_ts_pc(n,:,:) = ts_pc(:,:,stn_lat(n),stn_lon(n));
%         stn_pr_pc(n,:,:) = pr_pc(:,:,stn_lat(n),stn_lon(n));
    end

    save([DIR_NAME,'/rnd',num2str(NUM_STNS),'_prox',num2str(m),'.mat'], ...
         'stn_lat','stn_lon','stn_ts');
     % 'stn_synts','stn_synpr','stn_ts_pc','stn_pr_pc' if necessary
end

% Writing README file

fid = fopen([DIR_NAME,'/README.txt'], 'w+');

fprintf(fid,'%% This file contains the parameters that were used with station_select.m ');
fprintf(fid,['when generating the files in this folder(',DIR_NAME,').\n\n']);
fprintf(fid,'Boundaries of the box where the stations are contained.\n');
fprintf(fid,['S_lat = ',num2str(S_lat),'\n']);
fprintf(fid,['N_lat = ',num2str(N_lat),'\n']);
fprintf(fid,['W_lon = ',num2str(W_lon),'\n']);
fprintf(fid,['E_lon = ',num2str(E_lon),'\n\n']);
fprintf(fid,'Selection from areas with absolute correlation over a certain threshold.\n');
fprintf(fid,['MIN_COR = ',num2str(MIN_COR),'\n']);
fprintf(fid,['COR_WDW = ',num2str(COR_WDW(1)),':',num2str(COR_WDW(end)),'\n\n']);
fprintf(fid,['Number of pseudoproxy stations: ',num2str(NUM_STNS),'\n\n']);
fprintf(fid,'Station selection conditions:\n');
fprintf(fid,'abs(corr_ts)>MIN_COR\n');
% fprintf(fid,'squeeze(abs(mean(ts_pc(1,33:end,:,:)-ts_pc(2,33:end,:,:),2))) > 0.3\n');
% fprintf(fid,'nonstat_tsmap > 50\n');
% fprintf(fid,'lsfrac > 0.75\n');
fprintf(fid,'Tropical regions have not been included\n');

fclose(fid);



% %% Viewing avaliable Data on a Map
% 
% % Nonstationarities Map
% load nonstat_map.mat
% subplot(2,1,1);
% pcolor(lon,lat,nonstat_prmap);
% plotworld;
% colorbar;
% caxis([0, 200]);
% colormap(flipud(hot));
% title('Num of nonstationary stations for prec')
% subplot(2,1,2);
% pcolor(lon,lat,nonstat_tsmap);
% plotworld;
% caxis([0, 200]);
% colorbar
% colormap(flipud(hot));
% title('Num of nonstationary stations for temp')

% % EOF Plots
% load runcorr_eofs.mat
% load runcorr.mat
% NUM_OF_EOFS = 10;
% eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
% eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
% figure;
% for n=1:2
%     subplot(2,2,n)
%     contourf(lon,lat,squeeze(eof_ts_fm(n,:,:)));
%     plotworld;
%     caxis([-0.025,0.025])
%     colormap(redblue(13))
%     caxis([-0.025,0.025])
%     title(['EOF ',num2str(n),' of Temp running correlations, Explained Variance: ',num2str(expvar_ts(n)),'%']);
%     subplot(2,2,n+2)
%     plot(PC_ts(n,:))
% end
% figure;
% for n=1:2
%     subplot(2,2,n)
%     pcolor(lon,lat,squeeze(eof_pr_fm(n,:,:)));
%     plotworld;
%     colormap(redblue(13))
%     title(['EOF ',num2str(n),' of Precip running correlations, Explained Variance: ',num2str(expvar_pr(n)),'%']);
%     subplot(2,2,n+2)
%     plot(PC_pr(n,:))
% end
% %% Examining Selected Stations
% 
% % Showing Location of Stations
% RUN_NUM = 1;
% NUM_STNS = 3;
% DIR_NAME = 'Pseudoproxies_highlat_chstn' ;
% % load(['Pseudoproxies_nonpolar/rnd10_prox',num2str(RUN_NUM),'.mat'])
% load([DIR_NAME,'/',num2str(NUM_STNS),'stns_1000prox.mat']);
% axis([0 360 -90 90]);
% clf
% plotworld
% hold on
% scatter(lon(stn_lon(:)),lat(stn_lat(:)),'ro','filled')
% hold off
% %% Correlation Time Series
% 
% for n=1:NUM_STNS
%     subplot(ceil(NUM_STNS/2),2,n)
% %     plot SOMETHING
%     xlabel('End Year of Correlation Window')
%     ylabel('Correlation Coefficients')
%     title('');
% end

