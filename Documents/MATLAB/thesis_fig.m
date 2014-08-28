% This script will have all the figures to plot to the thesis

%% Setup
ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
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

%% Figure 1-1
figure;
window = 31; NUM_YRS = 499;
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

% Temp
bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_ts(:));
sorted_ts_runcorr = ts_runcorr(:,sorted_corr_ind); % This will need to be checked
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; ts_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    ts_runcorr_quan(:,m) = quantile(reshape(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% The plotting part

clf; axes; hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(HA([4]),'LineStyle','-','LineWidth',3);
set(text(0,0.30,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.07,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',30,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['Running Correlations (',num2str(window),' yr windows)'],'FontSize',14);
xlabel('Correlations over 499 yr period','FontSize',14);
title(['Correlation percentiles for modeled temperature, rcor=',num2str(window),'yrs'],'FontSize',14)
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.03 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -1:0.2:1, ...
    'XTick'       , -1:0.2:1, ...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1], ...
    'FontSize'    , 14   ...
    );
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm


%% Figure 1-2
figure;

contourf(lon,lat,corr_ts); axis equal;
shading flat
title('Correlation Coefficients for Temperature','FontSize',14,'FontWeight','bold');
plotworld;
colormap(b2r(-1,1));
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 14   ...
    );
colorbar;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 28 19]);

%% Figure 2-1
figure;

load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map31yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,1)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap31);
title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])
cbfreeze;
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map61yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,2)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap61);
title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])
cbfreeze;
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map91yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,3)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap91);
title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])
cbfreeze

%% Figure 2-2
figure;
for i=1:3
subplot(3,1,i)
window = i*30+1; NUM_YRS = 499;
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

% Temp
bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_ts(:));
sorted_ts_runcorr = ts_runcorr(:,sorted_corr_ind); % This will need to be checked
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; ts_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    ts_runcorr_quan(:,m) = quantile(reshape(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% The plotting part

hold on; Hnd = nan(1,7); cmap = hsv(8);
for n=7:-1:1
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color',cmap(n,:),'LineWidth',2);
end
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(HA([4]),'LineStyle','-','LineWidth',3);
ylabel(['Running Correlations (',num2str(window),' yr windows)']);
xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled temperature, rcor=',num2str(window),'yrs'])
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .01] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , [-1,-0.3,0,0.3,1], ...
    'XTick'       , [-1,-0.3,0,0.3,1], ...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1] ...
    );

end

legend('99th Percentile','95th Percentile','75th Percentile','50th Percentile','25th Percentile','5th Percentile','1st Percentile');

%% Figure 3
figure;

load('DataFiles/500yrCalWdw_meth_stats.mat')
qEPC = quantile(all_stn_corr_EPC_RV,[.05 .5 .95], 2);
qMRV = quantile(all_stn_corr_MRV,[.05 .5 .95], 2);
qCPS = quantile(all_stn_corr_CPS_RV,[.05 .5 .95], 2);

clf; axes; hold on; grid on;
Hnd(1,1) = plot(squeeze(qMRV(:,1)),'--g');
set(Hnd(1,1),'Color','g','MarkerFaceColor','g');
Hnd(1,3) = plot(squeeze(qMRV(:,3)),'--g');
Hnd(2,1) = plot(squeeze(qEPC(:,1)),'--r');
Hnd(2,3) = plot(squeeze(qEPC(:,3)),'--r');
Hnd(3,1) = plot(squeeze(qCPS(:,1)),'--b');
Hnd(3,3) = plot(squeeze(qCPS(:,3)),'--b');
Hnd(1,2) = plot(squeeze(qMRV(:,2)),'g','LineWidth',2);
Hnd(2,2) = plot(squeeze(qEPC(:,2)),'r','LineWidth',2);
Hnd(3,2) = plot(squeeze(qCPS(:,2)),'b','LineWidth',2);
hold off;
set(Hnd(1,[1,3]),'Color','g','MarkerFaceColor','g');
set(Hnd(2,[1,3]),'Color','r','MarkerFaceColor','r');
set(Hnd(3,[1,3]),'Color','b','MarkerFaceColor','b');
ylabel('Percentile Correlations with Nino3.4 index','FontSize',14  ); 
xlabel('Number of Stations used in reconstruction','FontSize',14  );
legend([Hnd(1:3,2); Hnd(3,1); Hnd(3,3);],'MRV Median','EPC\_RV Median','CPS\_RV Median', ...
       '5^t^h Percentile','95^t^h Percentile','location','southeast'                );
% title('Reconstructions from a global selection of pseudoproxies, using all 499 years of data')
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'LineWidth'   , 1,  ...    
    'YLim'        , [0 1], ...
    'XLim'        , [0 70], ...
    'FontSize'    , 14 ...
    );

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm


%% Figure 4-7

for window = [31, 61, 91]


GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
end

% Plotting EPC
subplot(3,3,1+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
ylabel(['Correlation (',num2str(window),'yrs of data)'])
if window==31 title(['EPC\_RV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting CPS
subplot(3,3,2+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==91 xlabel('Number of Stations included in reconstruction'); end
if window==31 title(['CPS\_RV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
% Plotting MRV
subplot(3,3,3+(floor(window/30)-1)*3)
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on
if window==31 title(['MRV']); end
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 

end
suptitle([strrep(GROUP_NAME,'_','\_'),' - Ranges of Correlation percentiles'])
set(gcf, 'PaperPosition', [0 0 28 19]);
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);

%% Figure 8 - Probabilities of Non-stationary Stations
figure;
GROUP_NAME = 'glb_ts'; NUM_TRIALS = 1000; c=1;
for window = [31]
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; cmap = hsv(NUM_GROUPS); i=1; dataline = nan(NUM_GROUPS,42);
    for group_size=[3,5,10,30,50,70]
        
        all_nstat_yrs = nan([1,NUM_TRIALS]);
        all_nstat_nstns = nan([1,NUM_TRIALS]);
        corrs = nan([1,NUM_TRIALS]);

        DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
        load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
         'stn_lat','stn_lon','indice_pool');

        % Compare the lat and lon with nonstat map, save no. of nonstat, save
        % nonstat year nums, and no. of nonstat points in reconstruction Plot
        % numbers of non-stat stations/av years against Reconstruction skill. Also
        % find percentages of the runs that have the numbers of the non-stat
        % prox/av yrs

        nstat_yrs = nan(size(stn_lat));
        for m=1:group_size
            for tr=1:NUM_TRIALS
                nstat_yrs(tr,m) = nonstat_tsmap(stn_lat(tr,m),stn_lon(tr,m));
            end
        end

        nstat_avyrs = mean(nstat_yrs,2);
        nstat_nstns = sum(nstat_yrs > ceil(0.1*(NUM_YRS-window)),2)/group_size;
        dataline(i,:)=histline(nstat_nstns,[-0.025:0.025:1.0]+0.025/2);
        i=i+1;
    end

        
        for n=1:NUM_GROUPS
            s_dataline = smooth(dataline(n,:)); hold on;
            plot([0:0.025:1.025],s_dataline,'Color',cmap(n,:),'LineWidth',2);
        end
        hold off;
        grid minor
        xlim([0 1]); ylim([0 500]); xlabel(' Proportion of non-stationary stations ','FontSize',14);
        legend('Group Size of 3','Group Size of 5','Group Size of 10','Group Size of 30','Group Size of 50','Group Size of 70');
        % h = findobj(gca,'Type','Patch');
        % set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
        ylabel('Number of reconstructions (out of 1000)','FontSize',14)
%         title('Proportion of nstat stns per reconstruction','FontSize',14);

end

%% Appendix Figure 1


ts_file = 'DataFiles/ts_A1.nc';
lat = nc_varget(ts_file,'lat'); lon = nc_varget(ts_file,'lon');
for window=[31]
    
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 499;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:NUM_CAL_WDW-1
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
folderlist = dir(['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/']);
folderlist = folderlist([8 11 16 17]+2);% Remove the ones you dont want
for i=1:(length(folderlist))
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',folderlist(i).name];
    lat_ind = []; lon_ind=[]; num_prox = zeros(NUM_CAL_WDW,1);
    for c=1:NUM_CAL_WDW
        load([DIR_NAME, '/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/3stns_1000prox.mat']);
        [lat_ind_tmp,lon_ind_tmp]=ind2sub([90 144],indice_pool);
        lat_ind = [lat_ind; lat_ind_tmp]; lon_ind = [lon_ind; lon_ind_tmp];
        num_prox(c) = length(indice_pool);
    end
    subplot(4,1,rem(i-1,4)+1)
    values = hist3([lon_ind lat_ind],{1:length(lon) , 1:length(lat)})';
    pcolor(lon,lat,values); plotworld; colormap(flipud(hot(10))); colorbar;
    xlim([0 360]); ylim([-90 90]);
    title([strrep(DIR_NAME(58:end),'_','\_'),' - Density with all CalWdws, mean no.: ',num2str(mean(num_prox),'%.0f'),', std: ',num2str(std(num_prox),'%.0f')]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
                
end


end