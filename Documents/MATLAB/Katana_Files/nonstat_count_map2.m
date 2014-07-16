% This script wil plot the number of nonstationarities in an area of the
% map from existing data of synthetic running correlations.
% This requires the synthetic correlations produced by synth_runcorr.m, but
% transformed by synth_pointform.m in order to produce a map of nonstationarities
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

%% Fitting Distributions and Obtaining percentiles

nonstat_prmap=nan(length(lat),length(lon),'single');
nonstat_tsmap=nan(length(lat),length(lon),'single');
running_nonstat_prmap=nan(499,length(lat),length(lon),'single');
running_nonstat_tsmap=nan(499,length(lat),length(lon),'single');
nonstat_tsmaprecord=zeros(length(n34_ind),length(lat),length(lon),'single');
nonstat_prmaprecord=zeros(length(n34_ind),length(lat),length(lon),'single');
pr_pc = nan(2,length(n34_ind),length(lat),length(lon),'single');
ts_pc = nan(2,length(n34_ind),length(lat),length(lon),'single');
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']);
runcorr_wdw = window; % Running correlation window in the runcorr.mat file

% Limits of box
S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['../Data/Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            % This takes 0.3 seconds per point
            load(['../Data/Synth_pointform/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat']);
            spotz_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) ); % Fishers Z Score
            spotz_pr = 0.5*log( (1+spot_pr)./(1-spot_pr) );
            pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
            ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
            pr_pc_spotz = prctile(spotz_pr,[2.5,97.5]);
            ts_pc_spotz = prctile(spotz_ts,[2.5,97.5]);
            pr_pc(:,:,j,i) = prctile(spot_pr,[2.5,97.5]);
            ts_pc(:,:,j,i) = prctile(spot_ts,[2.5,97.5]);
            nonstat_prmap(j,i)=length( find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)')   );
            nonstat_tsmap(j,i)=length( find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)')   );
            nonstat_tsmaprecord(find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)'),j,i) = 1;
            nonstat_prmaprecord(find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)'),j,i) = 1;
                                        
        else
            disp(['Data at ',num2str(lon(i)),'E',num2str(lat(j)),'N does not exist']);
        end
    end
end

save(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'],'nonstat_prmap','nonstat_tsmap', ...
     'pr_pc','ts_pc',...
     'nonstat_tsmaprecord','nonstat_prmaprecord','window');
 
% %% Plotting Map
% subplot(2,1,1);
% contourf(lon,lat,nonstat_prmap,10);
% plotworld;
% colorbar;
% caxis([0, 200]);
% colormap(flipud(hot));
% title('Num of nonstationary stations for prec')
% subplot(2,1,2);
% contourf(lon,lat,nonstat_tsmap,10);
% plotworld;
% caxis([0, 200]);
% colorbar
% colormap(flipud(hot));
% title('Num of nonstationary stations for temp')
% 
% %% Adding Difference in Percentiles Overlay

% subplot(2,1,1)
% hold on;
% % [c,h]=contour(lon,lat,squeeze(max(abs(squeeze(pr_pc(1,:,:,:)-pr_pc(2,:,:,:))))),[0.1:0.1:0.9],'b');
% [c,h]=contour(lon,lat,squeeze(abs(mean(pr_pc(1,33:end,:,:)-pr_pc(2,33:end,:,:),2))),[0.1:0.2:0.7],'b');
% clabel(c,h,'Color','b')
% hold off;
% subplot(2,1,2)
% hold on
% % [c,h]=contour(lon,lat,squeeze(max(abs(squeeze(ts_pc(1,:,:,:)-ts_pc(2,:,:,:))))),[0.1:0.1:0.9],'b');
% [c,h]=contour(lon,lat,squeeze(abs(mean(ts_pc(1,33:end,:,:)-ts_pc(2,33:end,:,:),2))),[0.1:0.2:0.7],'b');
% clabel(c,h,'Color','b')
% hold off;
% % %% Removing insignificant teleconnections for Temperature
% % 
% % sig_chang         e=ones(size(ts_pc,3),size(ts_pc,4));
% % for i=1:size(ts_pc,3)
% %     for j=1:size(ts_pc,4)
% %         if abs(max(ts_pc(1,:,i,j)-ts_pc(2,:,i,j)))<0.1
% %             sig_change(i,j)=nan; % No significant change
% %         end
% %     end
% % end
% %% Examining Certain Stations plotting Z-scores
% 
% x_lon = 226; [~,x_ind]= min(abs(lon-x_lon));
% y_lat = 1; [~,y_ind]= min(abs(lat-y_lat));
% 
% load(['Synth_runcorr/',num2str(lon(x_ind)),'E',num2str(lat(y_ind)),'N_syncorr.mat']);
% load('DataFiles/runcorr.mat');     
% 
% spotz_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) );
% spotz_pr = 0.5*log( (1+spot_pr)./(1-spot_pr) );
% pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
% ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
% pr_pc_spotz = prctile(spotz_pr,[2.5,97.5]);
% ts_pc_spotz = prctile(spotz_ts,[2.5,97.5]);
% 
% % Plot of percentiles and GFDL runcorrs
% subplot(2,1,1)
% plot(pr_pc_spotz'); hold on;
% plot(squeeze(pr_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% title(['Percentiles of Z scores of Runing correlations of pr at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% subplot(2,1,2)
% plot(ts_pc_spotz'); hold on;
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% title(['Percentiles of Z scores of Runing correlations of ts at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% 
% 
% %% Plotting map of running nonstat
% clf;
% for i=1:(length(ts_runcorr(:,1,1))+1-run_wdw-runcorr_wdw)
%     pcolor(lon,lat,squeeze(running_nonstat_tsmap(i,:,:))); plotworld; colorbar; caxis([0, 30]);
%     title(['Running nonstationary years in 30 year windows (temperature)',num2str(i)])
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'Plots/running_nonstat_ts.avi','FPS',8)

% %% Plotting map of nonstat time series
% clf;
% for i=1:(length(ts_runcorr(:,1,1)))
%     pcolor(lon,lat,squeeze(nonstat_tsmaprecord(i,:,:))); plotworld; colorbar; caxis([0 1]);
%     colormap(flipud(gray)); title(['Nonstationary Areas (temperature)',num2str(i)])
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'Plots/nonstat_ts.avi','FPS',8)

% %% Plotting Different groups of nonstationarities
% 
% load DataFiles/nonstat_map.mat
% load DataFiles/runcorr.mat
% 
% corr_pr = zeros(size(apr,2),size(apr,3));
% for i=1:size(apr,2)
%     for j=1:size(apr,3)
%         corr_pr(i,j) = corr(n34_ind,apr(:,i,j));
%     end
% end
% 
% corr_ts = zeros(size(ats,2),size(ats,3));
% for i=1:size(ats,2)
%     for j=1:size(ats,3)
%         corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
%     end
% end
% clear i j
% 
% run_corr_low = zeros(size(corr_ts));
% mean_corr_low = zeros(size(corr_ts));
% for i=1:size(apr,2)
%     for j=1:size(apr,3)
%         if (~isempty(find(abs(ts_runcorr(:,i,j)) < 0.3)) & ...
%                                  nonstat_tsmap(i,j) > 50 & ...
%                                  abs(corr_ts(i,j)) > 0.3       )
%             run_corr_low(i,j) = 1;
%         elseif (~isempty(find(abs(ts_runcorr(:,i,j)) > 0.3)) & ...
%                                  nonstat_tsmap(i,j) > 50 & ...
%                                  abs(corr_ts(i,j)) < 0.3       )
%              mean_corr_low(i,j) = 1;
%         end
%     end
% end
% 
% [is js] = ind2sub(size(corr_ts),find(run_corr_low));
% [is2 js2] = ind2sub(size(corr_ts),find(mean_corr_low));
% plotworld; hold on; run_hnd=scatter(lon(js), lat(is),'r.'); mean_hnd=scatter(lon(js2),lat(is2),'b.');
% xlim([0,360]); ylim([-90,90]); hold off;
% legend([run_hnd, mean_hnd],'abs(runcorr) > 0.3','corr\_ts > 0.3')
