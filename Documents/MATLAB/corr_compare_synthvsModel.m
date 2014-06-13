% This script will look at the correlations between the synthetic data and
% the model data
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r', 'plotworld', and the folder DataFiles
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

window = 31; % The running window in years

% Loading Synthetic Data ~110s

nu_ts_100=zeros(100,499,90,144);
nu_pr_100=zeros(100,499,90,144);
for n=1:100
    load(['Synth_Data/run',num2str(n),'syn.mat'])
    nu_ts_100(n,:,:,:) = nu_ts;
    nu_pr_100(n,:,:,:) = nu_pr;
end

%% Calculating Correlations

corr_ats_nu = zeros(100,size(ats,2),size(ats,3));
corr_apr_nu = zeros(100,size(ats,2),size(ats,3));
for i=1:size(ats,2) % ~10 min
    for j=1:size(ats,3) 
        for n=1:100
            corr_ats_nu(n,i,j) = corr(squeeze(nu_ts_100(n,2:499,i,j))',ats(2:499,i,j));
            corr_apr_nu(n,i,j) = corr(squeeze(nu_pr_100(n,2:499,i,j))',apr(2:499,i,j));
        end
    end
end

corr_nu_nu_ts = nan(100,size(ats,2),size(ats,3));
corr_nu_nu_pr = nan(100,size(apr,2),size(apr,3));
for i=1:size(ats,2) % ~10 min
    for j=1:size(ats,3) 
        for n=2:100
            corr_nu_nu_ts(n,i,j) = corr(squeeze(nu_ts_100(1,2:499,i,j))',squeeze(nu_ts_100(n,2:499,i,j))');
            corr_nu_nu_pr(n,i,j) = corr(squeeze(nu_pr_100(1,2:499,i,j))',squeeze(nu_pr_100(n,2:499,i,j))');
        end
    end
end

%% Plotting the correlations
figure;
subplot(3,1,1)
pcolor(lon,lat,squeeze(mean(corr_apr_nu,1)))
colorbar
colormap(b2r(-1,1));
plotworld
title('Correlation of Model with Synthetic Prec anomalies')
subplot(3,1,2)
pcolor(lon,lat,squeeze(mean(corr_nu_nu_pr(2:100,:,:),1)))
colorbar
colormap(b2r(-1,1));
plotworld
title('Correlation of Synthetic with Synthetic Prec anomalies')
subplot(3,1,3)
pcolor(lon,lat,squeeze(mean(corr_apr_nu,1)-mean(corr_nu_nu_pr(2:100,:,:),1)))
colorbar
colormap(b2r(-0.05,0.05));
plotworld
title('Subtraction of Middle plot from the top plot')

% Variance
figure;
subplot(3,1,1)
pcolor(lon,lat,squeeze(var(corr_apr_nu)))
colorbar
% colormap(b2r(0,5e-3));
caxis([0,5e-3]);
plotworld
title('Variance of Correlation of Model with Synthetic pr anomalies')
subplot(3,1,2)
pcolor(lon,lat,squeeze(var(corr_nu_nu_pr(2:100,:,:))))
colorbar
caxis([0,5e-3]);
% colormap(b2r(0,5e-3));
plotworld
title('Variance of Correlation of Synthetic with Synthetic pr anomalies')
subplot(3,1,3)
pcolor(lon,lat,squeeze(var(corr_apr_nu)-var(corr_nu_nu_pr(2:100,:,:))))
colorbar
% colormap(b2r(-1e-5,1e-5));
caxis([-1e-3,1e-3]);
plotworld
title('Subtraction of Middle plot from the top plot')