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
RV_WDW = [15:(499-14)]; c=1;

GROUP_NAME = 'glb_ts';
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];                              %%%%%%%
load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/10stns_1000prox'])
