% This script will use existing running correlations between synthetic temp/prec
% and the Nino3.4 index (from synth_corr), modifying the format so each
% grid point can be analysed for percentiles and nonstationarities

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

window = 31; % The running window in years

% Loading 
load('runcorr.mat');
pr_runcorrz = pr_runcorr;
pr_runcorrz = 0.5*log( (1+pr_runcorrz)./(1-pr_runcorrz) );
ts_runcorrz = ts_runcorr;
ts_runcorrz = 0.5*log( (1+ts_runcorrz)./(1-ts_runcorrz) );

%% Changing data organisation to spatial-point form
tic;
% Limits of box
S_lat = -90; N_lat = -60; W_lon = 150; E_lon = 180;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

spot_ts = zeros(1000,499);
spot_pr = zeros(1000,499);
for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['Synth_corr_',num2str(window),'yrwin_polar/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            load(['Synth_corr_',num2str(window),'yrwin_polar/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'])
            if size(spot_pr,1)<1000
                recalc=true;
            else 
                recalc=false;
            end
        else
            recalc=true;
        end
        if recalc==true
            for n=1:1000 % This will take about 6 minutes
                load(['Synth_corr_',num2str(window),'yrwin_polar/run',num2str(n),'syncorr.mat']);
                spot_ts(n,:) = ts_synruncorr(:,j,i);
                spot_pr(n,:) = pr_synruncorr(:,j,i);
            end
            save(['Synth_corr_',num2str(window),'yrwin_polar/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],...
                'spot_ts','spot_pr');
        end
        toc;
    end
end
toc;


%% Analysis of a single point
x_lon = 226; [~,x_ind]= min(abs(lon-x_lon));
y_lat = 3; [~,y_ind]= min(abs(lat-y_lat));

load(['Synth_corr_',num2str(window),'yrwin_0-120E/',num2str(lon(x_ind)),'E',num2str(lat(y_ind)),'N_syncorr.mat']);

% Translating to Fisher Z scores

spotz_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) );
spotz_pr = 0.5*log( (1+spot_pr)./(1-spot_pr) );

% Fitting Distributions and Obtaining percentiles of CORR

% for n=33:499
%     [pr_muhat pr_sighat] = normfit(spotz_pr(:,n));
%     [ts_muhat ts_sighat] = normfit(spotz_ts(:,n));
%     if kstest2(spotz_pr(:,n),normrnd(pr_muhat,pr_sighat,1,10000))
%         disp(['Warning: Normal Series fit (pr) is inappropriate at n=',num2str(n)]);
%     endcd
%     if kstest2(spotz_ts(:,n),normrnd(ts_muhat,ts_sighat,1,10000))
%         disp(['Warning: Normal Series fit (ts) is inappropriate at n=',num2str(n)]);
%     end    
% end

subplot(2,1,1)
pr_pc_spot = prctile(spot_pr,[2.5,97.5]);
ts_pc_spot = prctile(spot_ts,[2.5,97.5]);
plot(pr_pc_spot'); hold on;
plot(squeeze(pr_runcorr(:,y_ind,x_ind)),'k','LineWidth',3);
hold off;
title(['No. of Nonstationary Prec Stations at ', num2str(x_lon),'E, ',num2str(y_lat),'N is ',...
    num2str(length(find(squeeze(pr_runcorr(:,y_ind,x_ind))<pr_pc_spot(1,:)'|squeeze(pr_runcorr(:,y_ind,x_ind))>pr_pc_spot(2,:)')))])
ylabel('Correlations Coefs')
subplot(2,1,2)
plot(ts_pc_spot'); hold on;
plot(squeeze(ts_runcorr(:,y_ind,x_ind)),'k','LineWidth',3);
hold off;
title(['No. of Nonstationary Temp Stations at ', num2str(x_lon),'E, ',num2str(y_lat),'N is ',...
    num2str(length(find(squeeze(ts_runcorr(:,y_ind,x_ind))<ts_pc_spot(1,:)'|squeeze(ts_runcorr(:,y_ind,x_ind))>ts_pc_spot(2,:)')))])
ylabel('Correlations Coefs')
%% Plot of Z scores for several time series for precipitation

plot(squeeze(pr_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3)
hold on;
plot(squeeze(spotz_pr(1,:)),'r')
plot(squeeze(spotz_pr(2,:)),'y')
plot(squeeze(spotz_pr(3,:)),'b')
plot(squeeze(spotz_pr(4,:)),'g')
plot(squeeze(spotz_pr(5,:)),'m')
hold off
title(['Fisher Z scores of Runing correlations of GFDL data vs synthetic data at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
ylabel('Fisher Z Correlation Coefficients');
xlabel('Year after the final year of the running correlation period');
legend('GFDL','Synth1','Synth2','Synth3','Synth4','Synth5');

%% Plotting Area with data
lat = nc_varget(ts_file,'lat');
lon = nc_varget(ts_file,'lon');
load Synth_corr_31yrwin/run1000syncorr.mat
pcolor(lon,lat,squeeze(ts_synruncorr(33,:,:)))
title('Correlation Coefficients of last synthetic run, first correlation period');
shading flat
plotworld