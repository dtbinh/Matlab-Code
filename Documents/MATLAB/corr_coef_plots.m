% Script calculates the correlation coefficients of the temperature or
% rainfall against Nino3.4 SSTs
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


%% Correlations

corr_pr = zeros(size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        corr_pr(i,j) = corr(n34_ind,apr(:,i,j));
    end
end

corr_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
    end
end

%% Plotting
% Precip Correlation plot (whole year)
figure(1);
contourf(lon,lat,corr_pr);
shading flat
title('Correlation Coefficients for Precipitation','FontSize',14,'FontWeight','bold');
plotworld;
colormap(b2r(-1,1));
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360], ...
    'FontSize'    , 14   ...
    );
colorbar;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 28 19]);
saveas(gcf,'Plots/corr(n34,apr).jpg')
% Temp Correlation Plot (whole year)
figure(2);
contourf(lon,lat,corr_ts);
shading flat
title('Correlation Coefficients for Temperature','FontSize',14,'FontWeight','bold');
plotworld;
colormap(b2r(-1,1));
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 14   ...
    );
colorbar;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 28 19]);
saveas(gcf,'Plots/corr(n34,ats).jpg')