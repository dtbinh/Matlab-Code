% This script will calculate running correlations between synthetic temp/prec
% and the Nino3.4 index, and saves them
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

window = 61; % The running window in years

%% Calculating running correlations from SYNTHETIC data
tic;
% Limits of box to calculate corr coefs
S_lat = -90; N_lat = 90; 
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
% [~,W_bound]= min(abs(lon-W_lon));
% [~,E_bound]= min(abs(lon-E_lon));

for n=beginning:1000

load(['/srv/scratch/z3372730/Synth_Data/run',num2str(n),'syn.mat'])

% USE IF OLD FILES EXIST ALREADY
if  ~exist(['Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'],'file')
    ts_synruncorr=NaN(size(nu_ts),'single');
    pr_synruncorr=NaN(size(nu_pr),'single');
else
    load(['Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'])
end

% Running Correlation of Temperature
for i=S_bound:N_bound
    for j=W_bound:E_bound
        ts_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_ts(:,i,j)),n34_ind],window,2);
        % Note that this running correlation places the value after the window
    end
end

% Running Correlation of Precipitation
for i=S_bound:N_bound
    for j=W_bound:E_bound
        pr_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_pr(:,i,j)),n34_ind],window,2);
    end
end

save(['Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'],'ts_synruncorr','pr_synruncorr','window');
toc;
end

% %% Gathering Data for PDFs
% 
% x_lon = 160; [~,x_ind]= min(abs(lon-x_lon));
% y_lat = -20; [~,y_ind]= min(abs(lat-y_lat));
% 
% load(['Synth_corr_',num2str(window),'yrwin/run',num2str(1000),'syncorr.mat']);
% if sum(~isnan(pr_synruncorr(:,y_ind,x_ind)),1)==0
%     error('Selected spot has no data');
% end
% 
% spot_ts=zeros(1000,499); % 1000 runs, 499 years of corr data
% spot_pr=zeros(1000,499); % 1000 runs, 499 years of corr data
% 
% for n=1:1000 % This will take about 2.5 minutes
%     load(['Synth_corr_',num2str(window),'yrwin/run',num2str(n),'syncorr.mat']);
%     spot_ts(n,:) = squeeze(ts_synruncorr(:,y_ind,x_ind));
%     spot_pr(n,:) = squeeze(pr_synruncorr(:,y_ind,x_ind));
% end
% 
% % Translating to Fisher Z scores
% 
% spot_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) );
% spot_pr = 0.5*log( (1+spot_pr)./(1-spot_pr) );
% pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
% ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
% 
% %% Fitting Distributions and Obtaining percentiles
% 
% for n=33:499
%     [pr_muhat pr_sighat] = normfit(spot_pr(:,n));
%     if kstest2(spot_pr(:,n),normrnd(pr_muhat,pr_sighat,1,10000))
%         disp(['Warning: Normal Series fit is inappropriate at n=',num2str(n)]);
%     end
% end
% % Percentiles
% pr_pc_spot = prctile(spot_pr,[2.5,97.5]);
% ts_pc_spot = prctile(spot_ts,[2.5,97.5]);
% plot(ts_pc_spot'); hold on;
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% 
% length(find(squeeze(pr_runcorrz(:,y_ind,x_ind))<pr_pc_spot(1,:)'|squeeze(pr_runcorrz(:,y_ind,x_ind))>pr_pc_spot(2,:)'))
% length(find(squeeze(ts_runcorrz(:,y_ind,x_ind))<ts_pc_spot(1,:)'|squeeze(ts_runcorrz(:,y_ind,x_ind))>ts_pc_spot(2,:)'))
% %% Plot of Z scores for several time series
% 
% plot(squeeze(pr_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3)
% hold on;
% plot(squeeze(spot_pr(1,:)),'r')
% plot(squeeze(spot_pr(2,:)),'y')
% plot(squeeze(spot_pr(3,:)),'b')
% plot(squeeze(spot_pr(4,:)),'g')
% plot(squeeze(spot_pr(5,:)),'m')
% hold off
% title(['Fisher Z scores of Runing correlations of GFDL data vs synthetic data at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% ylabel('Fisher Z Correlation Coefficients');
% xlabel('Year after the final year of the running correlation period');
% legend('GFDL','Synth1','Synth2','Synth3','Synth4','Synth5');
% 
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3)
% hold on;
% plot(squeeze(spot_ts(1,:)),'r')
% plot(squeeze(spot_ts(2,:)),'y')
% plot(squeeze(spot_ts(3,:)),'b')
% plot(squeeze(spot_ts(4,:)),'g')
% plot(squeeze(spot_ts(5,:)),'m')
% hold off
% title(['Fisher Z scores of Runing correlations of GFDL data vs synthetic data at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% ylabel('Fisher Z Correlation Coefficients');
% xlabel('Year after the final year of the running correlation period');
% legend('GFDL','Synth1','Synth2','Synth3','Synth4','Synth5');
% 
% %% Plotting Area with data
% 
% load Synth_corr_31yrwin/run1000syncorr.mat
% pcolor(lon,lat,squeeze(ts_synruncorr(33,:,:)))
% title('Correlation Coefficients of last synthetic run, first correlation period');
% shading flat
% plotworld
