% This matlab script will produce synthetic running correlations for
% precipitation and temperature by using the same method as specified in
% {Gallant et al, 2013, 'Nonstationary Australasian Teleconnections and
%   Implications for Paleoclimate Reconstructions', Journal of Climate,
%   vol. 26 pp. 8827-8849} as equation (1): MODIFIED - no CONSTANT TERM

% nu(t) = a1*c(t) + sigma_nu*sqrt(1-r^2)*[eta_nu(t) + B*eta_nu(t-1)]

% where: nu(t) is the synthetic precipitation/temperature series
%        a0 is the first regression coefficient (also the mean of the data)
%        a1 is the second regression coefficient
%        c(t) is the Nino3.4 Index
%        sigma_nu is the standard deviation of the actual precip/temp series
%        r is the correlation between Nino3.4 index and precip/temp
%        eta_nu is random Gaussian noise
%        B is the autocorrelation of the climate variable at lag 1
%        [eta_nu(t) + B*eta_nu(t-1)] is red noise

% Note: This script requires the mexcdf package to be installed, and the
% following files: plotworld.m, coast_v2.mat, b2r.m to be in the current
% directory and in DataFiles

%% Setup
ts_file = 'DataFiles/ts_A1.nc';

lat = nc_varget(ts_file,'lat');
lon = nc_varget(ts_file,'lon');
time = nc_varget(ts_file,'time'); % Assumes both files use the same time
ts = nc_varget(ts_file,'ts')-273.15; % To Celsius
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nE] = min(abs(lon-240));
[~,nW] = min(abs(lon-190));
load('DataFiles/ella_ts.mat');
ats=ella_ts;
n34_ind = mean(mean(ats(:,nS:nN,nW:nE),3),2);

window = 31;

%% Calculating Running Correlations
% if ~exist('../DataFiles/ella_runcorr31yrwdw.mat','file')
% Limits of box to calculate corr coefs
S_lat = -15; N_lat = 15; W_lon = 130; E_lon = 300;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

% % Running Correlation of Temperature
% ts_runcorr=nan(size(ats));
% for i=S_bound:N_bound
%     for j=W_bound:E_bound
% 	ts_runcorr(:,i,j)=movingCorrelation([squeeze(ats(:,i,j)),n34_ind],window,2);
% 	% Note that this running correlation places the value after the window
%     end
% end
% 
% save(['DataFiles/ella_runcorr',num2str(window),'yrwdw.mat'],'ts_runcorr');
% display('ella_runcorr saved');
% 
% end

%% Regressions, Correlations, Autocorrelations and Standard Deviations
% reg_ats = zeros(size(ats,2),size(ats,3));
% for i=1:size(ats,2)
%     for j=1:size(ats,3)
%         reg_ats(i,j) = regress(ats(:,i,j),n34_ind);
%     end
% end
% 
% corr_ts = zeros(size(ats,2),size(ats,3));
% for i=1:size(ats,2)
%     for j=1:size(ats,3)
%         corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
%     end
% end
% 
% atcorr_ts = zeros(2,size(ats,2),size(ats,3));
% for i=1:size(ats,2)
%     for j=1:size(ats,3)
%         atcorr_ts(:,i,j) = autocorr(ats(:,i,j),1);
%     end
% end
% 
% atcorr_ts = squeeze(atcorr_ts(2,:,:));
% 
% 
% sigma_ts = zeros(size(ats,2),size(ats,3));
% for i=1:size(ats,2)
%     for j=1:size(ats,3)
%         sigma_ts(i,j) = std(ats(:,i,j));
%     end
% end
% assert(S_lat == -15);
% %% Calculating the Synthetic Series
% 
% ts_series = zeros(1000,499,90,144,'single');
% ts_synruncorr = nan(499,90,144);
% mkdir('../Data/ella_Synth_Data')
% mkdir(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/'])
% tic;
% for n=1:1000
%     eta_nu = randn(length(n34_ind),1);
%     nu_ts=zeros(length(n34_ind),size(ats,2),size(ats,3),'single');
%     nu_ts(1,:,:)=NaN;
% 
%     for i=S_bound:N_bound
%         for j=W_bound:E_bound
%             nu_ts(2:end,i,j) = reg_ats(i,j)*n34_ind(2:end) + ...
%                 sigma_ts(i,j)*sqrt(1.0-corr_ts(i,j)^2) * ...
%                 (eta_nu(2:end) + atcorr_ts(i,j)*eta_nu(1:end-1));
%         end
%     end
%     
%    save(['../Data/ella_Synth_Data/run',num2str(n),'syn.mat'],'nu_ts','eta_nu')
%    
%     % Running Correlation of Temperature
%     for i=S_bound:N_bound
%         for j=W_bound:E_bound
%             ts_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_ts(:,i,j)),n34_ind],window,2);
%             % Note that this running correlation places the value after the window
%         end
%     end
%     save(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'],'ts_synruncorr','window');
%     ts_series(n,:,:,:) = ts_synruncorr;
%     toc;
% end
% display('ella_synth_data and runcorr made');
% assert(S_lat == -15);

ts_series = zeros(1000,499,90,144,'single');
for n=1:1000
    load(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat']);
    ts_series(n,:,:,:) = ts_synruncorr;
end
display('Stuff reloaded')
    
mkdir(['../Data/ella_Synth_pointform/',num2str(window),'yrWindow/']);
for i=W_bound:E_bound
    for j=S_bound:N_bound
        spot_ts = ts_series(:,:,j,i);
        save(['../Data/ella_Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],...
            'spot_ts','window');
    end
end

display('ella_synthpointform complete')

%% Fitting Distributions and Obtaining percentiles

nonstat_tsmap=nan(length(lat),length(lon),'single');
running_nonstat_tsmap=nan(499,length(lat),length(lon),'single');
ts_pc = nan(2,length(n34_ind),length(lat),length(lon),'single');
load(['DataFiles/ella_runcorr',num2str(window),'yrwdw.mat']);
runcorr_wdw = window; % Running correlation window in the runcorr.mat file

for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['../Data/ella_Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            % This takes 0.3 seconds per point
            load(['../Data/ella_Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat']);
            spotz_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) ); % Fishers Z Score
            ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
            
            ts_pc_spotz = prctile(spotz_ts,[2.5,97.5]);
            
            ts_pc(:,:,j,i) = prctile(spot_ts,[2.5,97.5]);
            
            nonstat_tsmap(j,i)=length( find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)')   );            
                                        
        else
            disp(['Data at ',num2str(lon(i)),'E',num2str(lat(j)),'N does not exist']);
        end
    end
end

save(['DataFiles/ella_nonstat_map',num2str(window),'yrwdw.mat'],'nonstat_tsmap', ...
     'ts_pc',...
     'window');
     
display('ella_nonstat_map complete')
