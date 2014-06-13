% Script calculates the correlation coefficients of the temperature or
% rainfall against Nino3.4 SSTs
% This script needs mexcdf to be installed in matlab to run, and the
% function 'b2r'

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

%% Correlations

corr_pr = zeros(size(pra,2),size(pra,3));
for i=1:size(pra,2)
    for j=1:size(pra,3)
        corr_pr(i,j) = corr(n34_ind,pra(:,i,j));
    end
end

corr_ts = zeros(size(tsa,2),size(tsa,3));
for i=1:size(tsa,2)
    for j=1:size(tsa,3)
        corr_ts(i,j) = corr(n34_ind,tsa(:,i,j));
    end
end

corr_pr_mon = zeros(12,size(pra,2),size(pra,3));
for mon=1:12
    for i=1:size(pra,2)
        for j=1:size(pra,3)
            corr_pr_mon(mon,i,j) = corr(n34_ind(mon:12:end),pra(mon:12:end,i,j));
        end
    end
end

corr_ts_mon = zeros(12,size(tsa,2),size(tsa,3));
for mon=1:12
    for i=1:size(tsa,2)
        for j=1:size(tsa,3)
            corr_ts_mon(mon,i,j) = corr(n34_ind(mon:12:end),tsa(mon:12:end,i,j));
        end
    end
end

%% Plotting
% Precip Correlation plot (whole year)

pcolor(lon,lat,corr_pr)
shading flat
colorbar;
title('Correlation Coefficients of Nino3.4 and Anomalous Precipitation');
hold on;
load('coast_v2.mat');
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
hold off
colormap(b2r(-1,1));

% Precip Correlation plot (monthly)
clf
for mon=1:12
    subplot(4,3,mon);
    pcolor(lon,lat,squeeze(corr_pr_mon(mon,:,:)));
    shading flat
    colormap(jet);
    axis([0,360,-90,90]);
    caxis([-1,1]);
    title(['Month ',num2str(mon)]);
    hold on;
    load('coast_v2.mat');
    plot(long_e,lati_e,'k');
    plot(long_w,lati_w,'k');
    hold off
    colormap(b2r(-1,1));
end

axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
c=colorbar ('FontSize',15);
colormap(b2r(-1,1));
caxis([-1,1]);
ylabel(c,'Corr Coef')
suptitle('Corr coefs of Anomalous Rainfall on the Nino3.4 Index');

% Temp Correlation Plot (whole year)
pcolor(lon,lat,corr_ts);
shading flat;
colorbar;
title('Correlation Coefficients of Nino3.4 and Anomalous Temperature');
hold on;
load('coast_v2.mat');
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
hold off
colormap(b2r(-1,1));

% Temp Correlation Plot (Monthly)
clf
for mon=1:12
    subplot(4,3,mon);
    pcolor(lon,lat,squeeze(corr_ts_mon(mon,:,:)));
    shading flat
    colormap(jet);
    axis([0,360,-90,90]);
    caxis([-1,1]);
    title(['Month ',num2str(mon)]);
    hold on;
    load('coast_v2.mat');
    plot(long_e,lati_e,'k');
    plot(long_w,lati_w,'k');
    hold off
    colormap(b2r(-1,1));
end

axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
c=colorbar ('FontSize',15);
colormap(b2r(-1,1));
caxis([-1,1]);
ylabel(c,'Corr Coef')
suptitle('Corr coefs of Anomalous Temperature on the Nino3.4 Index');