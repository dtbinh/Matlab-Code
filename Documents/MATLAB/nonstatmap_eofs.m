% This script will compute the first few EOFs of the nonstationarities map
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

% Optional Loading

load('DataFiles/nonstat_map.mat')

%% Formatting for EOFs

NUM_OF_EOFS = 5;
% Weight according to latitude ? Nah

% Probably dont need to do it for our purposes

pr_runcorr_fm = reshape(pr_runcorr(32:end,:,:),size(pr_runcorr(32:end,:,:),1),size(pr_runcorr,2)*size(pr_runcorr,3));
ts_runcorr_fm = reshape(ts_runcorr(32:end,:,:),size(ts_runcorr(32:end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
tic;
[eof_pr,PC_pr,expvar_pr] = caleof(pr_runcorr_fm, NUM_OF_EOFS, 1);
[eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 1);
toc;

% save('runcorr_eofs.mat','eof_ts','eof_pr','PC_pr','PC_ts','expvar_pr','expvar_ts');
eof_pr_fm = reshape(eof_pr,NUM_OF_EOFS,size(pr_runcorr,2),size(pr_runcorr,3));
eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));

%% EOF Plots
figure;
for n=1:
    subplot(2,1,n)
    contourf(lon,lat,squeeze(eof_ts_fm(n,:,:)));
    plotworld;
    caxis([-0.025,0.025])
    colormap(redblue(13))
    caxis([-0.025,0.025])
    title(['EOF ',num2str(n),' of Temp running correlations, Explained Variance: ',num2str(expvar_ts(n)),'%']);
    subplot(2,3,n+3)
    plot(PC_ts(n,:))
end
figure;
for n=1:3
    subplot(2,3,n)
    pcolor(lon,lat,squeeze(eof_pr_fm(n,:,:)));
    plotworld;
    colorbar
    colormap(redblue(13))
    title(['EOF ',num2str(n),' of Precip running correlations, Explained Variance: ',num2str(expvar_pr(n)),'%']);
    subplot(2,3,n+3)
    plot(PC_pr(n,:))
end
