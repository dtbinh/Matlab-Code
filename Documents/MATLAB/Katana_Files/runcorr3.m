% This script calculates the running correlations of Temp/Precip against
% Nino 3.4 SSTs
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r',  'plotworld', and the folder DataFiles, as well as the function 'movingCorrelation'

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

%% Calculating Running Correlations

window = 91; % The running window in years
% Limits of box to calculate corr coefs
S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));


% Running Correlation of Temperature
ts_runcorr=zeros(size(ats));
for i=S_bound:N_bound
    for j=W_bound:E_bound
        ts_runcorr(:,i,j)=movingCorrelation([squeeze(ats(:,i,j)),n34_ind],window,2);
        % Note that this running correlation places the value after the window
    end
end

% Running Correlation of Precipitation
pr_runcorr=zeros(size(apr));
for i=S_bound:N_bound
    for j=W_bound:E_bound
        pr_runcorr(:,i,j)=movingCorrelation([squeeze(apr(:,i,j)),n34_ind],window,2);
    end
end
save(['DataFiles/runcorr',num2str(window),'yrwdw.mat'],'ts_runcorr','pr_runcorr');
pr_runcorr = pr_runcorr(window+1:end,:,:);
ts_runcorr = ts_runcorr(window+1:end,:,:);
pr_runcorr_mn = squeeze(mean(pr_runcorr,1));
ts_runcorr_mn = squeeze(mean(ts_runcorr,1));

% % Quantile Calculations
% 
% pr_quan = zeros(7,size(pr,2),size(pr,3));
% ts_quan = zeros(7,size(ts,2),size(ts,3));
% for i=1:size(pr,2)
%     for j=1:size(pr,3)
%         pr_quan(:,i,j) = quantile(squeeze(pr_runcorr(:,i,j)),[0,0.05,0.25,0.5,0.75,0.95,1]);
%         ts_quan(:,i,j) = quantile(squeeze(ts_runcorr(:,i,j)),[0,0.05,0.25,0.5,0.75,0.95,1]);
%     end
% end
% 
% %% Plot of Variability in Correlation
% % Rainfall
% figure;
% hold on;
% for i=1:size(pr_quan,1)
%     scatter(pr_runcorr_mn(:),pr_quan(i,:),1);
% end
% hold off;
% h=legend('0','0.05','0.25','0.5','0.75','0.95','1','Location','SouthEast');
% title(h,'Quantiles');
% axis equal
% axis([-1,1,-1,1]);
% grid on
% xlabel('Mean correlation for entire dataset');
% ylabel('30 year running correlation coefficients');
% title('Spread of correlation quantiles for Precipitation');
% 
% % Temperature
% figure;
% hold on;
% for i=1:size(ts_quan,1)
%     scatter(ts_runcorr_mn(:),ts_quan(i,:),1);
% end
% hold off;
% h=legend('0','0.05','0.25','0.5','0.75','0.95','1','Location','SouthEast');
% title(h,'Quantiles');
% axis equal
% axis([-1,1,-1,1]);
% grid on
% xlabel('Mean correlation for entire dataset');
% ylabel('30 year running correlation coefficients');
% title('Spread of correlation quantiles for Temperature');