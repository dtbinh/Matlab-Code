% This script calculates the running correlations of Temp/Precip against
% Nino 3.4 SSTs
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r',  'plotworld', and the folder DataFiles, as well as the function 'movingCorrelation'

%% Setup
load DataFiles/model_output.mat

%% Calculating Running Correlations

window = 31; % The running window in years
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
        tic;
        ts_runcorr(:,i,j)=movingCorrelation([squeeze(ats(:,i,j)),n34_ind],window,2);
        % Note that this running correlation places the value after the window
        toc;
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