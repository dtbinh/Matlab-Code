% This script will calculate the Atlantic Multidecadal Oscillation and the
% Indian Ocean Dipole index from Van Oldenburg et al. 2009 and Saji et al.
% 1999 respectively from the GFDL CM2.1 model data

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
time = nc_varget(ts_file,'time'); % Assumes both files use the same time
ts = nc_varget(ts_file,'ts')-273.15; % To Celsius
pr = nc_varget(pr_file,'pr');
load DataFiles/lsmask.mat

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

% Nino3.4
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nE] = min(abs(lon-240));
[~,nW] = min(abs(lon-190));
n34_ind = mean(mean(ats(:,nS:nN,nW:nE),3),2);

% Indian Ocean Dipole DMI
[~,iwN] = min(abs(lat-10));
[~,iwS] = min(abs(lat+10));
[~,iwE] = min(abs(lon-70));
[~,iwW] = min(abs(lon-50));
[~,ieN] = min(abs(lat-0));
[~,ieS] = min(abs(lat+10));
[~,ieE] = min(abs(lon-110));
[~,ieW] = min(abs(lon-90));
iod_ind = mean(mean(ats(:,iwS:iwN,iwW:iwE),3),2) - mean(mean(ats(:,ieS:ieN,ieW:ieE),3),2);

% AMO Index
[~,aN] = min(abs(lat-60));
[~,aS] = min(abs(lat-25));
[~,aE] = min(abs(lon-353));
[~,aW] = min(abs(lon-285));
ssta = nan(size(ats)); ssta(:,find(~lsmask)) = ats(:,find(~lsmask));

count = 0; runTotal = 0;
for i=aS:aN
    for j=aW:aE
        if ~isnan(ssta(1,i,j))
            runTotal = runTotal + ssta(:,i,j);
            count = count + 1;
        end
    end
end
amo_ind = runTotal./count;
clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

%% Correlations between running correlations and indices
window = 31;
load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

corr_ats_iod = nan(size(ats,2),size(ats,3));
corr_ats_amo = nan(size(ats,2),size(ats,3));
for i=1:length(lat)
    for j=1:length(lon)
        corr_ats_iod(i,j) = corr(ts_runcorr(window+1:end,i,j),iod_ind(window+1:end));
        corr_ats_amo(i,j) = corr(ts_runcorr(window+1:end,i,j),amo_ind(window+1:end));
    end
end  

CMAG = 0.4;
subplot(2,1,1)
pcolor(lon,lat,corr_ats_iod); plotworld;
title(['Correlation between iod-DMI and ts_runcorr(rcorwdw=',num2str(window),'r)']);
colorbar;
colormap(redblue(13))
caxis([-CMAG,CMAG])
subplot(2,1,2)
pcolor(lon,lat,corr_ats_amo); plotworld;
title(['Correlation between AMO and ts_runcorr(rcorwdw=',num2str(window),'r)']);
colorbar;
colormap(redblue(13))
caxis([-CMAG,CMAG])

% Time series comparison
amo = smooth(amo_ind,31);
iod = smooth(iod_ind,31);
iod = (iod-mean(iod))/std(iod);
amo = (amo-mean(amo))/std(amo);

for eof=1:NUM_OF_EOFS
    subplot(NUM_OF_EOFS,1,eof);
    hold on;
    plot(iod,'b')
    plot(amo,'g')
    plot(PC_ts_rcor(eof,:)/std(PC_ts_rcor(eof,:)),'k');
    grid on; hold off;
    iod_corr = corr(iod(16:end-16),squeeze(PC_ts_rcor(eof,:))');
    amo_corr = corr(amo(16:end-16),squeeze(PC_ts_rcor(eof,:))');
    title(['EOF',num2str(eof),' PC time series, iod corr: ',num2str(iod_corr,'%.2f'),', amo corr: ',num2str(amo_corr,'%.2f')]);
    ylim([-3,3]);
end
legend('IOD','AMO','PC\_time series')
corr(amo(16:end-16),iod(16:end-16));