% This script will compute the first few EOFs of the running correlations of the model data and show the
% variance explained
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

corr_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
    end
end

clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

window = 31; % The running window in years

% Optional Loading

load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])

%% Formatting for EOFs

NUM_OF_EOFS = 5;
% Weight according to latitude ? Nah

% Probably dont need to do it for our purposes

pr_runcorr_fm = reshape(pr_runcorr((window+1):end,:,:),size(pr_runcorr((window+1):end,:,:),1),size(pr_runcorr,2)*size(pr_runcorr,3));
ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
tic;
[eof_pr,PC_pr,expvar_pr] = caleof(pr_runcorr_fm, NUM_OF_EOFS, 2);
[eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 2);
toc;

% save('runcorr_eofs.mat','eof_ts_fm','eof_pr_fm','PC_pr','PC_ts','expvar_pr','expvar_ts');
eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));

%% EOF and PC Plots
figure;
for n=1:NUM_OF_EOFS
    subplot(NUM_OF_EOFS,2,2*n-1)
    contourf(lon,lat,squeeze(eof_ts_fm(n,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13)); colorbar;
    caxis([-0.025,0.025])
    title(['EOF ',num2str(n),' of Temp running correlations, Explained Variance: ',num2str(expvar_ts(n)),'%']);
    subplot(NUM_OF_EOFS,2,2*n)
    plot(PC_ts(n,:))
end
figure;
for n=1:NUM_OF_EOFS
    subplot(NUM_OF_EOFS,2,2*n-1)
    contourf(lon,lat,squeeze(eof_pr_fm(n,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13)); colorbar;
    caxis([-0.025,0.025])
    title(['EOF ',num2str(n),' of Prec running correlations, Explained Variance: ',num2str(expvar_pr(n)),'%']);
    subplot(NUM_OF_EOFS,2,2*n)
    plot(PC_pr(n,:))
end

%% EOF# Plots of all windows
THE_EOF = 10;
NUM_OF_EOFS = 10;
for window = [31,61,91]
    load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
    pr_runcorr_fm = reshape(pr_runcorr((window+1):end,:,:),size(pr_runcorr((window+1):end,:,:),1),size(pr_runcorr,2)*size(pr_runcorr,3));
    ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
    [eof_pr,PC_pr,expvar_pr] = caleof(pr_runcorr_fm, NUM_OF_EOFS, 2);
    [eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 2);
    eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
    eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
    subplot(3,2,2*floor((window/30))-1)
    contourf(lon,lat,squeeze(eof_ts_fm(1,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13)); colorbar;
    caxis([-0.025,0.025])
    title(['EOF ',num2str(THE_EOF),' of Temp, rcor=',num2str(window),'yr, Explained Variance: ',num2str(expvar_ts(THE_EOF)),'%']);
    subplot(3,2,2*floor((window/30)));
    contourf(lon,lat,squeeze(eof_pr_fm(1,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13)); colorbar;
    caxis([-0.025,0.025])
    title(['EOF ',num2str(THE_EOF),' of Prec, rcor=',num2str(window),'yr, Explained Variance: ',num2str(expvar_pr(THE_EOF)),'%']);
end


set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 28 19]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/eof',num2str(THE_EOF),'_rcor_pr&ts.jpg'])

% Explained variance plots
for window = [31,61,91]
    load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
    pr_runcorr_fm = reshape(pr_runcorr((window+1):end,:,:),size(pr_runcorr((window+1):end,:,:),1),size(pr_runcorr,2)*size(pr_runcorr,3));
    ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
    [eof_pr,PC_pr,expvar_pr] = caleof(pr_runcorr_fm, NUM_OF_EOFS, 2);
    [eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 2);
    eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
    eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
    subplot(3,2,2*floor((window/30))-1)
    plot(expvar_ts)
    grid on; ylim([0,40])
    title(['EOF Explained variance of Temp, rcor=',num2str(window),'yr']);
    subplot(3,2,2*floor((window/30)));
    plot(expvar_pr)
    title(['EOF Explained variance of Prec, rcor=',num2str(window),'yr']);
    grid on; ylim([0,40])
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/eof_expvar_rcor_pr&ts.jpg'])

%% EOF over correlation

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm

for n=1:3
    subplot(3,1,n)
    pcolor(lon,lat,corr_ts);
    plotworld;
    caxis([-1,1])
    colormap(b2r(-1,1)); colorbar;
    hold on;
    [cm, hand] = contour(lon,lat,squeeze(eof_ts_fm(n,:,:)),[-0.02,-0.01],'Color',[0.9 1.0 0.0],'LineWidth',2);
    clabel(cm,hand,'manual')
    [cm, hand] = contour(lon,lat,squeeze(eof_ts_fm(n,:,:)),[0.01,0.02],'Color',[0 1 0.2],'LineWidth',2);
    clabel(cm,hand,'manual')
    hold off;
    title(['Correlation of Temp plot under EOF',num2str(n) ,' contours. Explained Variance: ',num2str(expvar_ts(n)),'%']);
    saveas(gcf,['Plots/eof1-3_over_corr_ts.jpg'])
end

%% First five EOF plots of 31yrwdw runcorr

NUM_OF_EOFS = 5; window = 31;
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
pr_runcorr_fm = reshape(pr_runcorr((window+1):end,:,:),size(pr_runcorr((window+1):end,:,:),1),size(pr_runcorr,2)*size(pr_runcorr,3));
ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
[eof_pr,PC_pr,expvar_pr] = caleof(pr_runcorr_fm, NUM_OF_EOFS, 2);
[eof_ts,PC_ts_rcor,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 2);
eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
% EOF Maps
for THE_EOF=1:NUM_OF_EOFS
    subplot(ceil(NUM_OF_EOFS/2),2,THE_EOF)
    contourf(lon,lat,squeeze(eof_ts_fm(THE_EOF,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13));
    title(['EOF ',num2str(THE_EOF),' of Temp, rcor=',num2str(window),'yr, Explained Variance: ',num2str(expvar_ts(THE_EOF)),'%']);
end
% PC Timeseries
rgbmap = jet(5); 
for THE_EOF=1:NUM_OF_EOFS
    HA(THE_EOF) = plot(1:499-31,squeeze(PC_ts_rcor(THE_EOF,:)),'Color',rgbmap(THE_EOF,:)); hold on;
end
legend(HA,'PC1','PC2','PC3','PC4','PC5','location','southeast') hold off;
grid on
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 28*2 19*2]); %x_width=19cm y_width=28cm
% saveas(gcf,['Plots/eofs&pc_',num2str(window),'yr-rcor_ts.jpg'])

% Regression of above to global ts
reg_rcoreofs_ats = nan(NUM_OF_EOFS+1,size(ats,2),size(ats,3));
for i=1:length(lat)
    for j=1:length(lon)
        reg_rcoreofs_ats(:,i,j) = regress(squeeze(ats(1:499-window,i,j)), [PC_ts;ones(1,499-window)]');
    end
end

for eof=1:NUM_OF_EOFS+1
    subplot(ceil(NUM_OF_EOFS/2),2,eof)
    pcolor(lon,lat,squeeze(reg_rcoreofs_ats(eof,:,:)));
    plotworld;
    caxis([-0.07,0.07]);
    colormap(b2r(-0.07,0.07));
    title(['regress(ats,PC(caleof(rcor31yr))) regression coefficients, for EOF',num2str(eof)]);
end
title('Constant term from regression')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 28*2 19*2]); %x_width=19cm y_width=28cm
pause;
saveas(gcf,['Plots/reg(ats,eofs&pc_',num2str(window),'yr-rcor_ts).jpg'])