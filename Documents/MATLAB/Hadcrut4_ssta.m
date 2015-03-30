%% This script calculates the Nino3.4 teleconnection strengths on Hadcrut4 instrumental data

%% Date Conversions
START_YEAR = 1850;
START_MONTH = 1;
START_DAY = 1;
FILE_NAME = 'DataFiles/HadCRUT.4.3.0.0.median.nc';
DATE_OFFSET = datenum(START_YEAR,START_MONTH,START_DAY);
sst_datenum = nc_varget(FILE_NAME,'time')+DATE_OFFSET;


%% Nino 3.4 SSTs
ssta = nc_varget(FILE_NAME,'temperature_anomaly');
lat = nc_varget(FILE_NAME,'latitude');
lon = nc_varget(FILE_NAME,'longitude');
ssta_datenum = nc_varget(FILE_NAME,'time')+DATE_OFFSET;

% Nino 3.4 Indexes
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nW] = min(abs(lon-190));
[~,nE] = min(abs(lon-240));

% mnmean_sst=zeros(12,length(lat),length(lon));
% for n=1:12
%     mnmean_sst(n,:,:) = mean(mn_sst(n:12:end,:,:),1);
% end
% mn_ssta= zeros(size(mn_sst));
% for t=0:length(ssta_datenum)-1
%     mn_ssta(t+1,:,:) = mn_sst(t+1,:,:) - mnmean_sst(rem(t,12)+1,:,:);
% end

% Nino 3.4 Calculation

n34_ind = mean(mean(ssta(:,nN:nS,nW:nE),3),2);