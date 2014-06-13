% This script plots a quantile quantile plot of the model correlations over
% 500 years, vs the running correlation for the same period.

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

corr_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
    end
end
corr_pr = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_pr(i,j) = corr(n34_ind,apr(:,i,j));
    end
end

clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

load DataFiles/runcorr.mat
load DataFiles/nonstat_map.mat

%% Line Quantile Plots (newer version)


bin = -1.0:0.01:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_ts(:));
sorted_ts_runcorr = ts_runcorr(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1;
for m=1:length(bin_sizes)
    a = quantile(horzcat(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m)))),[0,0.05,0.25,0.5,0.75,0.95,1]);
    current_index = current_index + bin_sizes(m);
end




% pr_quan = zeros(7,size(apr,2),size(apr,3));
% ts_quan = zeros(7,size(ats,2),size(ats,3));
% for i=1:size(apr,2)
%     for j=1:size(apr,3)
%         pr_quan(:,i,j) = quantile(squeeze(pr_runcorr(:,i,j)),[0,0.05,0.25,0.5,0.75,0.95,1]);
%         ts_quan(:,i,j) = quantile(squeeze(ts_runcorr(:,i,j)),[0,0.05,0.25,0.5,0.75,0.95,1]);
%     end
% end
% 
% bin = -0.9:0.01:1.0;
% for n=1:7
%     [sorted_corr sorted_corr_ind] = sort(corr_ts(:));
%     sorted_ts_quan = ts_quan(n,sorted_corr_ind);
%     bin_sizes = histc(squeeze(sorted_corr),bin);
%     current_index = 1;
%     for m=1:length(bin_sizes)
%         mean_quan(n,m) = mean(sorted_ts_quan(current_index:(current_index-1+bin_sizes(m))));
%         current_index = current_index + bin_sizes(m);
%     end
% end

% The plotting part

axes; clf; hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,mean_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
hold off; grid on;
set(HA([1,7]),'LineStyle','-.','LineWidth',1);
set(text(0,-0.01,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
%% Plot of Variability in Correlation
% Rainfall
clf;
hold on;
for i=1:size(pr_quan,1)
    scatter(corr_pr(:),pr_quan(i,:),1,nonstat_prmap(:));
    colormap(flipud(hot))
    caxis([0 200])
end
hold off;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
axis equal
axis([-1,1,-1,1]);
grid on
xlabel('Mean correlation for entire dataset');
ylabel('30 year running correlation coefficients');
title('Spread of correlation quantiles for Precipitation');
h=colorbar;
title(h,'No. of Nonstationary Years')

% Temperature
figure;
clf;
hold on;
for i=1:size(ts_quan,1)
    scatter(corr_ts(:),ts_quan(i,:),1,nonstat_tsmap(:));
    colormap(flipud(hot))
    caxis([0 200])
end
hold off;
% h=legend('0','0.05','0.25','0.5','0.75','0.95','1','Location','SouthEast');
% title(h,'Percentiles');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
axis equal
axis([-1,1,-1,1]);
grid on
xlabel('Mean correlation for entire dataset');
ylabel('30 year running correlation coefficients');
title('Spread of correlation quantiles for Temperature');
h=colorbar;
title(h,'No. of Nonstationary Years')