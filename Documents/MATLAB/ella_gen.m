% This matlab script will produce synthetic running correlations for
% precipitation and temperature by using the same method as specified in
% {Gallant et al, 2013, 'Nonstationary Australasian Teleconnections and
%   Implications for Paleoclimate Reconstructions', Journal of Climate,
%   vol. 26 pp. 8827-8849} as equation (1):

% nu(t) = a0 + a1*c(t) + sigma_nu*sqrt(1-r^2)*[eta_nu(t) + B*eta_nu(t-1)]

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
load('DataFiles/ella_ts.mat');
ats=ella_ts;
n34_ind = mean(mean(ats(:,nS:nN,nW:nE),3),2);

window = 31;

%% Regressions, Correlations, Autocorrelations and Standard Deviations
reg_ats = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        reg_ats(i,j) = regress(ats(:,i,j),n34_ind);
    end
end

reg_apr = zeros(size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        reg_apr(i,j) = regress(apr(:,i,j),n34_ind);
    end
end

corr_pr = zeros(size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        corr_pr(i,j) = corr(n34_ind,apr(:,i,j));
    end
end

corr_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
    end
end

atcorr_pr = zeros(2,size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        atcorr_pr(:,i,j) = autocorr(apr(:,i,j),1);
    end
end

atcorr_ts = zeros(2,size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        atcorr_ts(:,i,j) = autocorr(ats(:,i,j),1);
    end
end
atcorr_pr = squeeze(atcorr_pr(2,:,:));
atcorr_ts = squeeze(atcorr_ts(2,:,:));

sigma_pr = zeros(size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        sigma_pr(i,j) = std(apr(:,i,j));
    end
end

sigma_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        sigma_ts(i,j) = std(ats(:,i,j));
    end
end

%% Calculating the Synthetic Series
mkdir('../Data/ella_Synth_Data')
for n=1:1000
    tic;
    eta_nu = randn(length(n34_ind),1);
    nu_ts=zeros(length(n34_ind),size(ats,2),size(ats,3),'single');
    nu_ts(1,:,:)=NaN;

    for i=1:size(ats,2)
        for j=1:size(ats,3)
            nu_ts(2:end,i,j) = ats_anmn(i,j) + reg_ats(i,j)*n34_ind(2:end) + ...
                sigma_ts(i,j)*sqrt(1.0-corr_ts(i,j)^2) * ...
                (eta_nu(2:end) + atcorr_ts(i,j)*eta_nu(1:end-1));
        end
    end

    nu_pr=zeros(length(n34_ind),size(apr,2),size(apr,3),'single');
    nu_pr(1,:,:)=NaN;

    for i=1:size(apr,2)
        for j=1:size(apr,3)
            nu_pr(2:end,i,j) = apr_anmn(i,j) + reg_apr(i,j)*n34_ind(2:end) + ...
                sigma_pr(i,j)*sqrt(1.0-corr_pr(i,j)^2) * ...
                (eta_nu(2:end) + atcorr_pr(i,j)*eta_nu(1:end-1));
        end
    end

   save(['../Data/ella_Synth_Data/run',num2str(n),'syn.mat'],'nu_ts','nu_pr','eta_nu')
   toc;
end

%% Calculating running correlations from SYNTHETIC data
tic;
% Limits of box to calculate corr coefs
S_lat = -90; N_lat = 90; 
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
% [~,W_bound]= min(abs(lon-W_lon));
% [~,E_bound]= min(abs(lon-E_lon));
mkdir(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/'])
for n=1:1000

load(['/srv/scratch/z3372730/Data/ella_Synth_Data/run',num2str(n),'syn.mat'])

% USE IF OLD FILES EXIST ALREADY
if  ~exist(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'],'file')
    ts_synruncorr=NaN(size(nu_ts),'single');
    pr_synruncorr=NaN(size(nu_pr),'single');
else
    load(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'])
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

save(['../Data/ella_Synth_runcorr/',num2str(window),'yrWindow/run',num2str(n),'syncorr.mat'],'ts_synruncorr','pr_synruncorr','window');
	
toc;
end
%% Power Spectral Density

% rednoise=(eta_nu(2:end) + atcorr_ts(i,j)*eta_nu(1:end-1));
% Fs=length(rednoise);
% t = 0:1/Fs:1-1/Fs;
% 
% N = length(rednoise);
% xdft = fft(rednoise);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)).*abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(rednoise):Fs/2;
% plot(freq,10*log10(psdx)); grid on;
% title('Periodogram Using FFT');
% xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');