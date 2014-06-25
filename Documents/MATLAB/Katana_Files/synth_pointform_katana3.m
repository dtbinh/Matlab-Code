% This script will use existing running correlations between synthetic temp/prec
% and the Nino3.4 index (from synth_corr), modifying the format so each
% grid point can be analysed for percentiles and nonstationarities
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r', 'plotworld', and the folder DataFiles

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';

lat = nc_varget(ts_file,'lat');
lon = nc_varget(ts_file,'lon');

window = 91; % The running window in years

%% Changing data organisation to spatial-point form

% Limits of box
S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

ts_series = zeros(1000,499,90,144,'single');
pr_series = zeros(1000,499,90,144,'single');
tic;

for n=1:1000 % This will take about 6 minutes ??????
    for m=1:24
        load(['../Data/',num2str(m),'Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat']);
        ts_series(n,:,:,((m-1)*6+1):(m*6)) = ts_synruncorr;
        pr_series(n,:,:,((m-1)*6+1):(m*6)) = pr_synruncorr;
        toc;
    end
end
mkdir(['Data/Synth_pointform/',num2str(window),'yrWindow/']);
for i=W_bound:E_bound
    for j=S_bound:N_bound
        spot_ts = ts_series(:,:,j,i);
        spot_pr = pr_series(:,:,j,i);
        save(['Data/Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],...
            'spot_ts','spot_pr','window');
    end
end