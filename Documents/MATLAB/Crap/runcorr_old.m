% This script calculates the running correlations of Temp/Precip against
% Nino 3.4 SSTs
% This script needs mexcdf to be installed in matlab to run, and the
% function 'b2r', as well as the function 'movingCorrelation'

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

ts_monmn = zeros(12,length(lat),length(lon)); % This will be 500yr mean
pr_monmn = zeros(12,length(lat),length(lon));
for month=1:12
    ts_monmn(month,:,:) = mean(ts(month:12:end,:,:),1);
    pr_monmn(month,:,:) = mean(pr(month:12:end,:,:),1);
end

tsa = zeros(size(ts));
pra = zeros(size(pr));
for t=0:length(time)-1
    tsa(t+1,:,:) = ts(t+1,:,:) - ts_monmn(rem(t,12)+1,:,:);
    pra(t+1,:,:) = pr(t+1,:,:) - pr_monmn(rem(t,12)+1,:,:);
end
n34_ind = mean(mean(tsa(:,nS:nN,nW:nE),3),2);

%% Calculating Running Correlations

window = 31*12; % The running window in months
% Limits of box to calculate corr coefs
S_lat = -25; N_lat = -20; W_lon = 140; E_lon = 145;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));
% Warning, one moving correlation takes about 0.5-1 second


% Running Correlation of Temperature
ts_runcorr=zeros(size(tsa));
for i=S_bound:N_bound
    for j=W_bound:E_bound
        ts_runcorr(:,i,j)=movingCorrelation([squeeze(tsa(:,i,j)),n34_ind],window,2);
    end
end

% Running Correlation of Precipitation
pr_runcorr=zeros(size(pra));
for i=S_bound:N_bound
    for j=W_bound:E_bound
        pr_runcorr(:,i,j)=movingCorrelation([squeeze(pra(:,i,j)),n34_ind],window,2);
    end
end

%% Plotting

% Location of Running Correlations
figure(1);
load('coast_v2.mat');
hold on;
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
axis([0,360,-90,90]);
plot([W_lon W_lon E_lon E_lon W_lon],[S_lat N_lat N_lat S_lat S_lat],'r');
hold off
title('Location of Correlations area');


% Plots of all running Correlations (tsa)
figure(2);
clf
hold on;
for i=S_bound:N_bound
    for j=W_bound:E_bound
        plot(ts_runcorr(:,i,j));
    end
end
hold off
title(['Running Corr of tsa and N34with ',num2str(window/12), ' year window, and boundaries: ', ...
    num2str(W_lon),' ',num2str(E_lon),' ',num2str(S_lat),' ',num2str(N_lat)]);

% Plot of a single running correlation
figure(3);
clf
corr_lat_ind = S_bound;
corr_lon_ind = E_bound;
plot(ts_runcorr(:,corr_lat_ind, corr_lon_ind));
title(['Running Corr of tsa and N34 with ',num2str(window/12), ' year window at: ', ...
    num2str(lat(corr_lat_ind)),' Lat, ',num2str(lon(corr_lon_ind)),' Lon']);
        
% Plots of all running Correlations (pra)
figure(4);
clf
hold on;
for i=S_bound:N_bound
    for j=W_bound:E_bound
        plot(pr_runcorr(:,i,j));
    end
end
hold off
title(['Running Corr of pra and N34 with ',num2str(window/12), ' year window, and boundaries: ', ...
    num2str(W_lon),' ',num2str(E_lon),' ',num2str(S_lat),' ',num2str(N_lat)]);

% Plot of a single running correlation
figure(5);
clf

corr_lat_ind = S_bound;
corr_lon_ind = E_bound;
plot(pr_runcorr(:,corr_lat_ind, corr_lon_ind));
title(['Running Corr of pra and N34 with ',num2str(window/12), ' year window at: ', ...
    num2str(lat(corr_lat_ind)),' Lat, ',num2str(lon(corr_lon_ind)),' Lon']);

% To remove NaNs A(isnan(A(:)))=[]