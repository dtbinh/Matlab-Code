% This script will find and plot teleconnections to the Nino3.4 index
% generated from the model data, supposedly
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
    
%% Teleconnections
telcon_ts = zeros(size(tsa,2),size(tsa,3));
for i=1:size(tsa,2)
    for j=1:size(tsa,3)
        telcon_ts(i,j) = regress(tsa(:,i,j),n34_ind);
    end
end

telcon_ts_mon = zeros(12,size(tsa,2),size(tsa,3));
for mon=1:12
    for i=1:size(tsa,2)
        for j=1:size(tsa,3)
            telcon_ts_mon(mon,i,j) = regress(tsa(mon:12:end,i,j),n34_ind(mon:12:end));
        end
    end
end

% Temp Reg Plot
for mon=1:12
    subplot(4,3,mon);
    set(gca, 'color', [0 0 0]);
    hold on;
    pcolor(lon,lat,squeeze(telcon_ts_mon(mon,:,:)));
    hold off;
    shading flat
    colormap(jet);
    axis([0,360,-90,90]);
    caxis([-2,2]);
    title(['Month ',num2str(mon)]);
    hold on;
    load('coast_v2.mat');
    plot(long_e,lati_e,'k');
    plot(long_w,lati_w,'k');
    hold off
end

axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
c=colorbar ('FontSize',15);
colormap(b2r(-2,2));
caxis([-2,2]);
ylabel(c,'Regression Coefficients')
suptitle('Regression coefficients of Anomalous Temperature on the Nino3.4 Index');

% Temp Reg plot (annual)
pcolor(lon,lat,telcon_ts);
shading flat;
colorbar;
colormap(b2r(-1,1.5));
title('Regression Coefficients of Anomalous Temperature on the N34 Index');
hold on;
load('coast_v2.mat');
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
hold off

% Precip Reg
telcon_pr = zeros(size(pr,2),size(pr,3));
for i=1:size(pr,2)
    for j=1:size(pr,3)
        telcon_pr(i,j) = regress(pra(:,i,j),n34_ind);
    end
end

telcon_pr_mon = zeros(12,size(pr,2),size(pr,3));
for mon=1:12
    for i=1:size(pr,2)
        for j=1:size(pr,3)
            telcon_pr_mon(mon,i,j) = regress(pra(mon:12:end,i,j),n34_ind(mon:12:end));
        end
    end
end

% Precip Reg Plot
% for mon=1:12
%     pcolor(lon,lat,squeeze(telcon_pr_mon(mon,:,:)));
%     shading flat;
%     colorbar;
%     caxis([min(telcon_pr_mon(:)),max(telcon_pr_mon(:))]);
%     colormap(flipud(b2r(-0.01,0.01)));
%     pause
% end

for mon=1:12
    subplot(4,3,mon);
    pcolor(lon,lat,squeeze(zscore(telcon_pr_mon(mon,:,:))));
    shading flat
    axis([0,360,-90,90]);
    caxis([-4,4]);
    title(['Month ',num2str(mon)]);
    hold on;
    load('coast_v2.mat');
    plot(long_e,lati_e,'k');
    plot(long_w,lati_w,'k');
    hold off
end
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
c=colorbar ('FontSize',15);
colormap(b2r(-4,4));
caxis([-4,4]);
ylabel(c,'Regression Coefficients')
suptitle('Zscore of Regression coefficients of Anomalous Precipitation on the Nino3.4 Index');

% Precip Reg Plot (annual)
pcolor(lon,lat,telcon_pr);
shading flat;
colorbar;
title('Regression Coefficients of Anomolous Rainfall on N34 Index');
hold on;
load('coast_v2.mat');
plot(long_e,lati_e,'k');
plot(long_w,lati_w,'k');
hold off
colormap(b2r(-2*10^-5,2*10^-5))

%% Plotting

% Temperature Plot
for t=1:6000
    figure(1);
    pcolor(lon,lat,squeeze(ts(t,:,:)))
    shading flat
    h=colorbar;
    caxis([-30,max(max(max((ts))))]);
    title(h,'Temp {\circ}C')
    title(num2str(n34_ind(t)));
    pause
end


% figure(2);
% pcolor(lon(nW:nE),lat(nS:nN),squeeze(tsa(1,nS:nN,nW:nE)))
% shading flat
% h=colorbar;
% title(h,'Temp {\circ}C')

figure(1);
plot(time/365,n34_ind)
title('Nino3.4 Index');
grid on

% Precipitation Plot
for t=1:6000
    figure(2);
    pcolor(squeeze(log_pr(t,:,:)).*lsmask);
    colormap(flipud(jet));
    shading flat
    h=colorbar;
    caxis([-30,-5]);
    title(h,'Log(pr)')
    title(['Day ',num2str(time(t))]);
    pause
    
end


