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
% window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];

%% Beggining of Loop

DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/glb_ts'] ; mkdir(DIR_NAME);
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
% Calibration windows set to being 10 overlapping windows over 499 years
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
% CAL_WDW = [1:50; 51:100; 101:150; 151:200; 201:250; 251:300; 301:350; 351:400; 401:450; 450:499];

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

for c=1:size(CAL_WDW,1)

    corr_ts = nan(size(ats,2),size(ats,3),'single');
    for i=S_bound:N_bound
        for j=W_bound:E_bound
            corr_ts(i,j) = corr(n34_ind(CAL_WDW(c,:)),ats(CAL_WDW(c,:),i,j));
        end
    end
    
    mkdir([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end))])

for NUM_STNS = 3:STN_MAX
    
    stn_ts = nan(NUM_TRIALS,NUM_STNS, NUM_YRS,'single');
    stn_pr = nan(NUM_TRIALS,NUM_STNS, NUM_YRS,'single');
    stn_lat = nan(NUM_TRIALS,NUM_STNS);
    stn_lon = nan(NUM_TRIALS,NUM_STNS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conditions - Replace 1 with desired conditions

    indice_pool = find(abs(corr_ts)>MIN_COR & ...
                       1 & ...
                       1 & ...
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
    
    save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
             'stn_lat','stn_lon','stn_ts','indice_pool','corr_ts','window');

end

% Writing README file

fid = fopen([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/README.txt'], 'w+');

fprintf(fid,'%% This file contains the parameters that were used with corr_vs_NumStns_gen.m ');
fprintf(fid,['when generating the files in this folder(',DIR_NAME,').\n\n']);
fprintf(fid,'Boundaries of the box where the stations are contained.\n');
fprintf(fid,['S_lat = ',num2str(S_lat),'\n']);
fprintf(fid,['N_lat = ',num2str(N_lat),'\n']);
fprintf(fid,['W_lon = ',num2str(W_lon),'\n']);
fprintf(fid,['E_lon = ',num2str(E_lon),'\n\n']);
fprintf(fid,'Selection from areas with absolute correlation over a certain threshold.\n');
fprintf(fid,['MIN_COR = ',num2str(MIN_COR),'\n']);
fprintf(fid,['CAL_WDW = ',num2str(CAL_WDW(c,1)),':',num2str(CAL_WDW(c,end)),'\n\n']);
fprintf(fid,'Data is also using temperature only\n');
% fprintf(fid,'Data is also using precipitation only\n');
fprintf(fid,'This file was produced using the UNSW Katana Computational Cluster.\n\n');
fprintf(fid,'Station selection conditions:\n');
fprintf(fid,'abs(corr_ts)>MIN_COR\n');
% fprintf(fid,'squeeze(abs(mean(ts_pc(1,33:end,:,:)-ts_pc(2,33:end,:,:),2))) > 0.3\n');
% fprintf(fid,'nonstat_tsmap > 50\n');
% fprintf(fid,'lsfrac > 0\n');
% fprintf(fid,'Tropical regions have not been included\n');
% fprintf(fid,'eof_ts_fm2 > 0.01\n');
fclose(fid);

end