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

NUM_YRS=499; NUM_TRIALS=1000; numstnstocompare=3:70;
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
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',35,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['r(',num2str(window),'yr)'],'FontSize',14);
xlabel('r(499yr)','FontSize',14);
title(['Percentiles'],'FontSize',14,'FontWeight','bold')
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
title('r(ts,Nino3.4)','FontSize',14,'FontWeight','bold');
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
load DataFiles/nstat_cmaps.mat
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map31yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,1)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap31);
num_nonstat = length(find(nonstat_tsmap>47));
title(['r(',num2str(window),'yr)'],'FontSize',14,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 14   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',14)
cbfreeze;
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map61yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,2)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap61);
num_nonstat = length(find(nonstat_tsmap>44));
title(['r(',num2str(window),'yr)'],'FontSize',14,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 14   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',14)
cbfreeze;
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map91yrwdw.mat
NUM_CONTOURS = 10;

subplot(3,1,3)
contourf(lon,lat,nonstat_tsmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(ncmap91);
num_nonstat = length(find(nonstat_tsmap>41));
title(['r(',num2str(window),'yr)'],'FontSize',14,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 14   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',14)
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
current_index = 1; ts_runcorr_quan = nan(5,length(bin_sizes));
for m=1:length(bin_sizes)
    ts_runcorr_quan(:,m) = quantile(reshape(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.5,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% Density Overlay
nonstat_map_ind = find(nonstat_tsmap > ceil(0.1*(NUM_YRS-window)) & 1 );
[a,b] = ind2sub(size(nonstat_tsmap),nonstat_map_ind);
corr_3d = permute(repmat(corr_ts,[1 1 499]),[3,1,2]);
runcr = nan(size(ts_runcorr)); corr_3d_fmt = nan(size(corr_3d));
for j=1:length(a)
    indices = find(nonstat_tsmaprecord(:,a(j),b(j)));
    corr_3d_fmt(indices,a(j),b(j)) = corr_3d(indices,a(j),b(j));
    runcr(indices,a(j),b(j)) = ts_runcorr(indices,a(j),b(j));
end
data1 = corr_3d_fmt;
data2 = runcr;
values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1});
cmap=hot(64);
colormap(flipud(cmap(24:end,:)))

% The plotting part

hold on; Hnd = nan(1,5); cmap = hsv(6);
for n=5:-1:1
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
grid on; axis equal; axis([-1 1 -1 1]);
set(HA([1,5]),'LineStyle','-','LineWidth',1); % hold off;
set(HA([3]),'LineStyle','-','LineStyle','-.');
set(HA([1,2]),'LineStyle','-','LineStyle','--');
ylabel(['r(',num2str(window),'yr)'],'FontSize'    , 14);
xlabel('r(499yr)','FontSize'    , 14);
% title(['Correlation percentiles for modeled temperature, rcor=',num2str(window),'yrs'])
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
    'XLim'        , [-1 1], ...
    'FontSize'    , 14 ...
    );

pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
colorbar
caxis([0,100])
end

h=legend('99th','95th','50th','5th','1st','Orientation','vertical');
v = get(h,'title');
set(v,'string','Percentiles','FontSize',14);

%% Figure 3
figure;

load('DataFiles/500yrCalWdw_meth_stats.mat')
qEPC = quantile(all_stn_corr_EPC_RV,[.05 .5 .95], 2);
qMRV = quantile(all_stn_corr_MRV,[.05 .5 .95], 2);
qCPS = quantile(all_stn_corr_CPS_RV,[.05 .5 .95], 2);
qRVM = quantile(all_stn_corr_RVM,[.05 .5 .95], 2);

clf; axes; hold on; grid on;
Hnd(1,1) = plot(squeeze(qMRV(:,1)),'--c');
set(Hnd(1,1),'Color','g','MarkerFaceColor','c');
Hnd(1,3) = plot(squeeze(qMRV(:,3)),'-.c');
Hnd(2,1) = plot(squeeze(qEPC(:,1)),'--r');
Hnd(2,3) = plot(squeeze(qEPC(:,3)),'-.r');
Hnd(3,1) = plot(squeeze(qCPS(:,1)),'--k');
Hnd(3,3) = plot(squeeze(qCPS(:,3)),'-.k');
Hnd(4,1) = plot(squeeze(qRVM(:,1)),'--m');
Hnd(4,3) = plot(squeeze(qRVM(:,3)),'-.m');
Hnd(1,2) = plot(squeeze(qMRV(:,2)),'c','LineWidth',2);
Hnd(2,2) = plot(squeeze(qEPC(:,2)),'r','LineWidth',2);
Hnd(3,2) = plot(squeeze(qCPS(:,2)),'k','LineWidth',2);
Hnd(4,2) = plot(squeeze(qRVM(:,2)),'m','LineWidth',2);
hold off;
set(Hnd(1,[1,3]),'Color','c','MarkerFaceColor','c');
set(Hnd(2,[1,3]),'Color','r','MarkerFaceColor','r');
set(Hnd(3,[1,3]),'Color','k','MarkerFaceColor','k');
set(Hnd(4,[1,3]),'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
set(Hnd(4,[2]),'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
ylabel('r(499yr)','FontSize',14  ); 
xlabel('Network Size','FontSize',14  );
legend([Hnd(1:4,2); Hnd(3,1); Hnd(3,3);],'MRV Median','EPC\_RV Median','CPS\_RV Median', 'RVM Median', ...
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

% % Statisical Difference in Methods
% subplot(3,1,1)
% ks_stat=nan(70,1); mw_stat=nan(70,1);
% for g=numstnstocompare
%     [~,ks_stat(g)]=kstest2(squeeze(all_stn_corr_MRV(g,:)),squeeze(all_stn_corr_EPC_RV(g,:)));
%     [mw_stat(g),~]=ranksum(all_stn_corr_MRV(g,:),all_stn_corr_EPC_RV(g,:));
% end
% semilogy(ks_stat,'b'); hold on; semilogy(mw_stat,'g'); hold off; ylim([10^-60,1]); xlim([3 70]); grid on;
% title('MRV vs EPC\_RV');
% set(gca,'YTick',[10^-60,10^-50,10^-40,10^-30,10^-20,10^-10,10^-2])
% 
% subplot(3,1,2)
% ks_stat=nan(70,1); mw_stat=nan(70,1);
% for g=numstnstocompare
%     [~,ks_stat(g)]=kstest2(squeeze(all_stn_corr_MRV(g,:)),squeeze(all_stn_corr_CPS_RV(g,:)));
%     [mw_stat(g),~]=ranksum(all_stn_corr_MRV(g,:),all_stn_corr_CPS_RV(g,:));
% end
% semilogy(ks_stat,'b'); hold on; semilogy(mw_stat,'g'); hold off; ylim([10^-60,1]); xlim([3 70]); grid on;
% title('MRV vs CPS\_RV');
% set(gca,'YTick',[10^-60,10^-50,10^-40,10^-30,10^-20,10^-10,10^-2])
% 
% subplot(3,1,3)
% ks_stat=nan(70,1); mw_stat=nan(70,1);
% for g=numstnstocompare
%     [~,ks_stat(g)]=kstest2(squeeze(all_stn_corr_CPS_RV(g,:)),squeeze(all_stn_corr_EPC_RV(g,:)));
%     [mw_stat(g),~]=ranksum(all_stn_corr_CPS_RV(g,:),all_stn_corr_EPC_RV(g,:));
% end
% semilogy(ks_stat,'b'); hold on; semilogy(mw_stat,'g'); hold off; ylim([10^-60,1]); xlim([3 70]); grid on;
% title('CPS\_RV vs EPC\_RV');
% set(gca,'YTick',[10^-60,10^-50,10^-40,10^-30,10^-20,10^-10,0.01])
% 
% % RVM and MRV performance - Figure 3 discrepancy with McGregor 2013 Fig4
% 
% NUM_STNS=3;
% scatter(all_stn_corr_MRV(NUM_STNS,:),all_stn_corr_RVM(NUM_STNS,:))
% ylim([-0.25 1]); xlim([-0.25 1]); grid on;
% xlabel('MRV'); ylabel('RVM');
% hold on; plot(0:0.1:1,0:0.1:1,'k'); hold off;


%% Figure 4-7
figure;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
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
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
    temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
end

% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r(',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle([strrep(GROUP_NAME,'_','\_')])
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);
% saveas('gcf','../../Dropbox/Literature/Writing/THESIS/AMS LaTeX Package v4.3.1/Figs/Fig4.pdf','output','pdf');

%% Figure 8 - Probabilities of Non-stationary Stations
figure;
GROUP_NAME = 'glb_ts'; NUM_TRIALS = 1000;
for window = [31]
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; i=1; NUM_CAL_WDW = 10; cmap = hsv(NUM_GROUPS); cmap(2,2)=0.8;
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
    end
    for group_size=[3,5,10,30,50,70]
        corrs = nan([1,NUM_TRIALS]);
        nstat_nstns = nan(NUM_CAL_WDW,NUM_TRIALS);
        for c=1:size(CAL_WDW,1)

            DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'stn_lat','stn_lon','indice_pool');

            nstat_yrs = nan(size(stn_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs(tr,m) = nonstat_tsmap(stn_lat(tr,m),stn_lon(tr,m));
                end
            end

            nstat_avyrs = mean(nstat_yrs,2);
            nstat_nstns(c,:) = sum(nstat_yrs > ceil(0.1*(NUM_YRS-window)),2)/group_size;
            
        end
        xbins = (0:1.0/group_size:1)-0.0000001; xbins(1)=0; xbins(end)=1;
        h=histc(nstat_nstns(:),xbins)/100;
        hold on; Hnd(i) = plot(xbins,h,'Color',cmap(i,:),'LineWidth',2);
        h(~h) = NaN;
        plot(xbins,h,'o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerSize',7);
        i=i+1;
    end
        
    hold off;
    grid minor
    xlim([0 1]); ylim([0 70]); xlabel(' Proportion of non-stationary stations ','FontSize',14);
    legend([Hnd],'Group Size of 3','Group Size of 5','Group Size of 10','Group Size of 30','Group Size of 50','Group Size of 70');
    % h = findobj(gca,'Type','Patch');
    % set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    ylabel('Percentage of reconstructions','FontSize',14)
%     title('Proportion of nstat stns per reconstruction','FontSize',14);

end

%% Figure 9 - Histogram Comparisons
figure;
window = 31; numstnstocompare=70;
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
s_Hnd = tight_subplot(3,4,[0.01 0.02],[0.1 0.05],[0.1 0.01]);

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
glb_temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
glb_temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
glb_temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
glb_temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
ntrop_temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
ntrop_temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
ntrop_temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
ntrop_temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    glb_temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    glb_temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    glb_temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
    glb_temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
end

GROUP_NAME = 'ntrop_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    ntrop_temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    ntrop_temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    ntrop_temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
    ntrop_temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
end

i=1; 
GROUPS = [5 20 50];
for group_size = GROUPS
% Plotting EPC
axes(s_Hnd(1+(i-1)*4))
cmap = hsv(size(CAL_WDW,1)+1); bins = 0:0.05:1; 
glb_Hnd = zeros(size(CAL_WDW,1),1); ntrop_Hnd = glb_Hnd;
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(glb_temp_corr_EPC_RV(group_size,c,:)),bins)/1000;
%     hold on; glb_Hnd(c) = plot(bins,h,'','Color',cmap(9,:),'LineWidth',2);
end
jbfill([bins],max(h,[],1),min(h,[],1),'b','k',[],0.5);
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(ntrop_temp_corr_EPC_RV(group_size,c,:)),bins)/1000;
%     hold on; ntrop_Hnd(c) = plot(bins,h,'-','Color',cmap(1,:),'LineWidth',1);
end
glb_all = squeeze(glb_temp_corr_EPC_RV(group_size,:,:));
ntrop_all = squeeze(ntrop_temp_corr_EPC_RV(group_size,:,:));
sig_diff = kstest2(glb_all(:),ntrop_all(:),0.0001)
ranksum(glb_all(:),ntrop_all(:),'alpha',0.01)
jbfill([bins],max(h,[],1),min(h,[],1),'y','k','add',0.5);

hold off; grid on; xlim([0 1]); ylim([0 0.7]);
ylabel(['Percentage of reconstructions (grp\_size=',num2str(group_size),')'])
% set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7],'XTick',[0:0.2:1]); 

% Plotting CPS
axes(s_Hnd(2+(i-1)*4))
glb_Hnd = zeros(size(CAL_WDW,1),1); ntrop_Hnd = glb_Hnd;
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(glb_temp_corr_CPS_RV(group_size,c,:)),bins)/1000;
end
glb_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'b','k',[],0.5);
hold off; grid on; xlim([0 1]); ylim([0 0.7]);
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(ntrop_temp_corr_CPS_RV(group_size,c,:)),bins)/1000;
end
glb_all = squeeze(glb_temp_corr_CPS_RV(group_size,:,:));
ntrop_all = squeeze(ntrop_temp_corr_CPS_RV(group_size,:,:));
sig_diff = kstest2(squeeze(glb_all(1,:)),squeeze(ntrop_all(1,:)),0.0001)
ranksum(glb_all(:),ntrop_all(:),'alpha',0.01)
ntrop_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'y','k','add',0.5);
% if group_size==50 xlabel('Correlation with Nino3.4','FontSize',14); end
% set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7],'XTick',[0:0.2:1]); 

% Plotting MRV
axes(s_Hnd(3+(i-1)*4))
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(glb_temp_corr_MRV(group_size,c,:)),bins)/1000;
end
glb_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'b','k',[],0.5);
% ks_matrix_bin = nan(10,10); ks_matrix_val = ks_matrix_bin;
% for i=1:10
%     for j=1:10
%         [ks_matrix_bin(i,j),ks_matrix_val(i,j)]=kstest2(squeeze(ntrop_all(i,:)),squeeze(ntrop_all(j,:)),0.0001);
%     end
% end

h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(ntrop_temp_corr_MRV(group_size,c,:)),bins)/1000;
end
glb_all = squeeze(glb_temp_corr_MRV(group_size,:,:));
ntrop_all = squeeze(ntrop_temp_corr_MRV(group_size,:,:));
sig_diff = kstest2(glb_all(:),ntrop_all(:),0.0001)
ranksum(glb_all(:),ntrop_all(:),'alpha',0.01)
ntrop_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'y','k','add',0.5);
hold off; grid on; xlim([0 1]); ylim([0 0.7]);
% set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7],'XTick',[0:0.2:1]); 

% Plotting RVM
axes(s_Hnd(4+(i-1)*4))
h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(glb_temp_corr_RVM(group_size,c,:)),bins)/1000;
end
glb_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'b','k',[],0.5);
% ks_matrix_bin = nan(10,10); ks_matrix_val = ks_matrix_bin;
% for i=1:10
%     for j=1:10
%         [ks_matrix_bin(i,j),ks_matrix_val(i,j)]=kstest2(squeeze(ntrop_all(i,:)),squeeze(ntrop_all(j,:)),0.0001);
%     end
% end

h=zeros(10,length(bins));
for c=1:size(CAL_WDW,1)
    h(c,:)=histc(squeeze(ntrop_temp_corr_RVM(group_size,c,:)),bins)/1000;
end
glb_all = squeeze(glb_temp_corr_RVM(group_size,:,:));
ntrop_all = squeeze(ntrop_temp_corr_RVM(group_size,:,:));
sig_diff = kstest2(glb_all(:),ntrop_all(:),0.0001)
ranksum(glb_all(:),ntrop_all(:),'alpha',0.01)
ntrop_Hnd=jbfill([bins],max(h,[],1),min(h,[],1),'y','k','add',0.5);
hold off; grid on; xlim([0 1]); ylim([0 0.7]);
% set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7],'XTick',[0:0.2:1]); 

i=i+1;
end; clear i;

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:0.7],'XTick', [0:0.2:1]); 
end

axes(s_Hnd(10)); xlabel('Correlation with Nino3.4','FontSize',14);
for i=1:3
    axes(s_Hnd(1+(i-1)*4));
    set(gca,'YTickLabel',[0:0.1:1]);
    ylabel(['Network Size: ',num2str(GROUPS(i))])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:0.2:1]);
end
% suptitle(['PDF range of glb and ntrop experiments'])
set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7]); 
legend([glb_Hnd, ntrop_Hnd],'glb','ntrop','Location','northwest')



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

% Proportion of tropical stations
load('/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/31yrWindow/ntrop_ts/CalWdw:1-31/3stns_1000prox.mat')

[~,S_bd]= min(abs(lat--10));
[~,N_bd]= min(abs(lat-10));
[~,W_bd]= min(abs(lon-100));
[~,E_bd]= min(abs(lon-300));
[lat_ind_tmp,lon_ind_tmp]=ind2sub([90 144],indice_pool);
trop_ind = find(lat_ind_tmp<N_bd & lat_ind_tmp>S_bd & lon_ind_tmp<E_bd & lon_ind_tmp>W_bd);

length(trop_ind)/length(indice_pool) % Proportion of tropical stations
% Location Check
plotworld; hold on; scatter(lon(lon_ind_tmp(trop_ind)),lat(lat_ind_tmp(trop_ind))); hold off;


%% Appendix Figure 2
figure;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
for window = [31, 61, 91]

GROUP_NAME = 'glb_pr'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
    temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
end

% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r(',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle([strrep(GROUP_NAME,'_','\_')])
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);

%% Appendix Figure 3
figure;
numstnstocompare=3:70;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
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
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_rmse_EPC_RV','all_stn_rmse_CPS_RV','all_stn_rmse_MRV','all_stn_rmse_RVM')
    temp_corr_EPC_RV(:,c,:) = all_stn_rmse_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_rmse_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_rmse_MRV;
    temp_corr_RVM(:,c,:) = all_stn_rmse_RVM;
end

% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on


% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.5]); grid on


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:0.5],'XTick', [0:20:70]); 
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.1:0.5]);
    ylabel(['RMSE(',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle([strrep(GROUP_NAME,'_','\_')])
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);
%% Appendix Figure 4
figure;
% window=31;
% NUM_OF_EOFS = 10;
% load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']); % This section will take 4600 seconds
% ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
% [eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 1);
% rc31_eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
% rc31_PC_ts = PC_ts;
% rc31_expvar = expvar_ts;
% window=61;
% load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']);
% ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
% [eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 1);
% rc61_eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
% rc61_expvar = expvar_ts;
% rc61_PC_ts = PC_ts;
% window=91;
% load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']);
% ts_runcorr_fm = reshape(ts_runcorr((window+1):end,:,:),size(ts_runcorr((window+1):end,:,:),1),size(ts_runcorr,2)*size(ts_runcorr,3));
% [eof_ts,PC_ts,expvar_ts] = caleof(ts_runcorr_fm, NUM_OF_EOFS, 1);
% rc91_eof_ts_fm = reshape(eof_ts,NUM_OF_EOFS,size(ts_runcorr,2),size(ts_runcorr,3));
% rc91_expvar = expvar_ts;
% rc91_PC_ts = PC_ts;
% save('DataFiles/all_rcor_ts_eofs.mat','rc91_eof_ts_fm','rc61_eof_ts_fm','rc31_eof_ts_fm','rc31_expvar','rc61_expvar','rc91_expvar', ...
%        'rc31_PC_ts','rc61_PC_ts','rc91_PC_ts');

load('DataFiles/all_rcor_ts_eofs.mat')
subplot(3,1,2)
pcolor(lon,lat,squeeze(rc31_eof_ts_fm(1,:,:))); plotworld; colormap(b2r(-0.03,0.03))
title('EOF1 (31yr window)'); colorbar;
% subplot(3,1,2)
% pcolor(lon,lat,squeeze(rc61_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% title('EOF 2 with window length 61 years')
% subplot(3,1,3)
% pcolor(lon,lat,squeeze(rc91_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% title('EOF 2 with window length 91 years')
% subplot(2,2,4)
subplot(3,1,1)
cmap=hsv(3);
plot(rc31_expvar,'','Color',cmap(1,:),'LineWidth',2); hold on;
plot(rc61_expvar,'','Color',cmap(2,:),'LineWidth',2);
plot(rc91_expvar,'','Color',cmap(3,:),'LineWidth',2); hold off;
ylim([0 50]); xlim([1,10]); legend('31 year window','61 year window','91 year window');
ylabel('Percentage %'); grid on;
xlabel('Number of EOF')
title('Explained Variance');



% Correlation between PC timeseries
% a=rc61_PC_ts(2,:)'; b=rc91_PC_ts(2,:)';
% a=rc61_eof_ts_fm(2,:)'; b=rc91_eof_ts_fm(2,:)';
% corr(a,b)
% cor = 0;
% for i=0:30
%     cor(i+1) = corr(a((i+1):400+i),b(1:400))
% end
% plot(cor); grid minor

% Plotting PC Time series
subplot(3,1,3)
a=rc31_PC_ts(1,:)'; b=rc61_PC_ts(1,:)'; c=rc91_PC_ts(1,:)';
plot(a,'','Color',cmap(1,:)); hold on;
plot(b,'','Color',cmap(2,:)); grid on;
plot(c,'','Color',cmap(3,:)); hold off
legend('31 year window','61 year window','91 year window','location','southeast')
xlabel('Year')
title('Principal Component 1')

%% Appendix Figure 5
figure;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
for window = [31, 61, 91]

GROUP_NAME = 'pneof_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_corr_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_corr_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_corr_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_corr_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_corr_MRV(:,c,:) = all_stn_corr_MRV;
    temp_corr_RVM(:,c,:) = all_stn_corr_RVM;
end

% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_corr_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,1]); grid on


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r(',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle([strrep(GROUP_NAME,'_','\_')])
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);