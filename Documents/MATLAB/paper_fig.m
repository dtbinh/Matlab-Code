% This script will have all the figures to plot to the thesis

%% Setup
tic;
load DataFiles/model_output.mat

VAR_WDW = 30; % Moving window for moving variance is 30 Years
window = 31; % The running window in years
n34_ind_RV = movingvar(n34_ind,VAR_WDW);
RV_WDW = [15:(499-14)];

corr_ts = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_ts(i,j) = corr(n34_ind,ats(:,i,j));
    end
end

NUM_YRS=499; NUM_TRIALS=1000; numstnstocompare=3:70;


%% Figure 1-1
close
figure;
subplot(2,2,1)
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

hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',3);
end
plot([-0.3,-0.3],[-1 1],'k','LineWidth',2); plot([0.3,0.3],[-1 1],'k','LineWidth',2)
plot([-1 1],[-0.3,-0.3],'k','LineWidth',2); plot([-1 1],[0.3,0.3],'k','LineWidth',2)
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',2);
set(HA([4]),'LineStyle','-','LineWidth',4);
set(text(0,0.30,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.07,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',35,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['r(',num2str(window),'yr)'],'FontSize',16);
xlabel('r(499yr)','FontSize',16);
title(['Percentiles'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.05 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -1:0.5:1, ...
    'XTick'       , -1:0.5:1, ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1], ...
    'FontSize'    , 16   ...
    );
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [19 19]);
set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm
text(-1+0.08,1-0.18,'a)','FontSize',22,'FontWeight','bold');

% Figure 1-2
subplot(2,2,2);

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
    'LineWidth'   , 2,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );
colorbar;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 19]);
set(gcf, 'PaperSize', [35 13]);
text(0+1,90-1,'b)','FontSize',22,'FontWeight','bold');

% print -painters -dpdf -r600 test.pdf
set(gcf, 'PaperSize', [35 13]);

% Figure 1-3


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
title(['r(',num2str(window),'yr)'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',18)
text(0+1,90-4,'a)','FontSize',22,'FontWeight','bold');
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
title(['r(',num2str(window),'yr)'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',18)
text(0+1,90-4,'b)','FontSize',22,'FontWeight','bold');
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
title(['r(',num2str(window),'yr)'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 1,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );
text(310,-40,['N=',num2str(num_nonstat)],'FontSize',18)
text(0+1,90-4,'c)','FontSize',22,'FontWeight','bold');
cbfreeze

%% Figure 2-2
for i=1:3
subplot(3,1,1*i)
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
plot([-0.3,-0.3],[-1 1],'k','LineWidth',2); plot([0.3,0.3],[-1 1],'k','LineWidth',2)
plot([-1 1],[-0.3,-0.3],'k','LineWidth',2); plot([-1 1],[0.3,0.3],'k','LineWidth',2)
grid on; axis equal; axis([-1 1 -1 1]);
set(HA([1,5]),'LineStyle','-','LineWidth',2); % hold off;
set(HA([3]),'LineStyle','-','LineStyle','-.');
set(HA([1,2]),'LineStyle','-','LineStyle','--');
ylabel(['r(',num2str(window),'yr)'],'FontSize'    , 16);
xlabel('r(499yr)','FontSize'    , 16);
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
    'LineWidth'   , 2,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1], ...
    'FontSize'    , 16 ...
    );

pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
caxis([0,100])
text(-1+0.01,1-0.01,'a)','FontSize',22,'FontWeight','bold');
end

colorbar('southoutside')
h=legend('99th','95th','50th','5th','1st','Orientation','vertical');
v = get(h,'title');
set(v,'string','Percentiles','FontSize',16);
% set(gcf, 'PaperPosition', [0 0 20 24]);
% set(gcf, 'PaperSize', [24 24]);
% saveas(gcf,'Plots/nonstatmap_perc.pdf','pdf')

%% Figure 3
figure;

load('DataFiles/500yrCalWdw_meth_stats.mat')
qEPC = quantile(all_stn_corr_EPC_RV,[.05 .5 .95], 2);
qMRV = quantile(all_stn_corr_MRV,[.05 .5 .95], 2);
qCPS = quantile(all_stn_corr_CPS_RV,[.05 .5 .95], 2);
qRVM = quantile(all_stn_corr_RVM,[.05 .5 .95], 2);

clf; axes; hold on; grid on;
Hnd(1,1) = plot(squeeze(qMRV(:,1)),'--c','LineWidth',2);
set(Hnd(1,1),'Color','g','MarkerFaceColor','c');
Hnd(1,3) = plot(squeeze(qMRV(:,3)),'-.c','LineWidth',2);
Hnd(2,1) = plot(squeeze(qEPC(:,1)),'--r','LineWidth',2);
Hnd(2,3) = plot(squeeze(qEPC(:,3)),'-.r','LineWidth',2);
Hnd(3,1) = plot(squeeze(qCPS(:,1)),'--k','LineWidth',2);
Hnd(3,3) = plot(squeeze(qCPS(:,3)),'-.k','LineWidth',2);
Hnd(4,1) = plot(squeeze(qRVM(:,1)),'--m','LineWidth',2);
Hnd(4,3) = plot(squeeze(qRVM(:,3)),'-.m','LineWidth',2);
Hnd(1,2) = plot(squeeze(qMRV(:,2)),'c','LineWidth',3);
Hnd(2,2) = plot(squeeze(qEPC(:,2)),'r','LineWidth',3);
Hnd(3,2) = plot(squeeze(qCPS(:,2)),'k','LineWidth',3);
Hnd(4,2) = plot(squeeze(qRVM(:,2)),'m','LineWidth',3);
hold off;
set(Hnd(1,[1,3]),'Color','c','MarkerFaceColor','c');
set(Hnd(2,[1,3]),'Color','r','MarkerFaceColor','r');
set(Hnd(3,[1,3]),'Color','k','MarkerFaceColor','k');
set(Hnd(4,[1,3]),'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
set(Hnd(4,[2]),'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
ylabel('r(499yr)','FontSize',16  ); 
xlabel('Network Size','FontSize',16  );
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
    'LineWidth'   , 2,  ...    
    'YLim'        , [0 1], ...
    'XLim'        , [0 70], ...
    'FontSize'    , 16 ...
    );

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 10]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [12 12]);

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


%% Figure 4-7 - Correlation as Skill Metric
clf
letters = 'abcdefghijkl';
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
in_range_499yr_pc_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
in_range_499yr_pc_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
in_range_499yr_pc_MRV = nan(max(numstnstocompare),length([31 61 91]));
in_range_499yr_pc_RVM = nan(max(numstnstocompare),length([31 61 91]));

for window = [31, 61, 91]

GROUP_NAME = 'glb_ts_nstat'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
% DIR_NAME = ['/home/nfs/z3372730/Documents/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

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

for n=numstnstocompare
    in_range_499yr_pc_EPC_RV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_corr_EPC_RV(n,:,:))>qEPC(n,1) & squeeze(temp_corr_EPC_RV(n,:,:))<qEPC(n,3) )))/(NUM_CAL_WDW*NUM_TRIALS);
    in_range_499yr_pc_CPS_RV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_corr_CPS_RV(n,:,:))>qCPS(n,1) & squeeze(temp_corr_CPS_RV(n,:,:))<qCPS(n,3) )))/(NUM_CAL_WDW*NUM_TRIALS);
    in_range_499yr_pc_MRV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_corr_MRV(n,:,:))>qMRV(n,1) & squeeze(temp_corr_MRV(n,:,:))<qMRV(n,3) )))/(NUM_CAL_WDW*NUM_TRIALS);
    in_range_499yr_pc_RVM(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_corr_RVM(n,:,:))>qRVM(n,1) & squeeze(temp_corr_RVM(n,:,:))<qRVM(n,3) )))/(NUM_CAL_WDW*NUM_TRIALS);
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
hold on; plot(qEPC,'k'); plot(squeeze(in_range_499yr_pc_EPC_RV(:,floor(window/30))),'k','LineWidth',2); hold off;
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
hold on; plot(qCPS,'k'); plot(squeeze(in_range_499yr_pc_CPS_RV(:,floor(window/30))),'k','LineWidth',2); hold off;
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
hold on; plot(qMRV,'k'); plot(squeeze(in_range_499yr_pc_MRV(:,floor(window/30))),'k','LineWidth',2); hold off;
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
hold on; plot(qRVM,'k'); plot(squeeze(in_range_499yr_pc_RVM(:,floor(window/30))),'k','LineWidth',2); hold off;
xlim([0,70]); ylim([0,1]); grid on

end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
    hold on; plot([0,70],[0.37 0.37],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
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

% suptitle([strrep(GROUP_NAME,'_','\_')])
suptitle('NSTAT_{glb\_ts}');
set(gcf, 'PaperPosition', [0 0 20 23]);
legendH = legend('5^t^h Percentile','95^t^h Percentile','Median','location','best','Orientation','horizontal');
set(legendH, 'FontSize',16,'Position',[0.372689213508714 0.0116864099741356 0.54468339307049 0.0419381107491857]);
saveas(gcf,['../../Dropbox/calWdwCorr_vs_NumStns_',GROUP_NAME,'.jpg']);

%% Figure 4-7 - Explained Variance as Skill Metric -Home
clf
letters = 'abcdefghijkl';

% Skilful Threshold
skilful_threshold = sqrt(0.5); % In correlation

s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
skilful_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
skilful_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
skilful_MRV = nan(max(numstnstocompare),length([31 61 91]));
skilful_RVM = nan(max(numstnstocompare),length([31 61 91]));

for window = [31, 61, 91]

GROUP_NAME = 'ntrop_ts_nstat'; % Change group name to get other figs
% DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
% DIR_NAME = ['/home/nfs/z3372730/Documents/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
temp_xvar_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_xvar_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_xvar_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_xvar_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_xvar_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_xvar_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_xvar_MRV(:,c,:) = all_stn_corr_MRV;
    temp_xvar_RVM(:,c,:) = all_stn_corr_RVM;
end

for n=numstnstocompare
    skilful_EPC_RV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_xvar_EPC_RV(n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
    skilful_CPS_RV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_xvar_CPS_RV(n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
    skilful_MRV(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_xvar_MRV(n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
    skilful_RVM(n,floor(window/30)) = sum(sum(sum(...
        squeeze(temp_xvar_RVM(n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
end
% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_xvar_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
hold on; plot(squeeze(skilful_EPC_RV(:,floor(window/30))),'k','LineWidth',2); hold off;
xlim([0,70]); ylim([0,1]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_xvar_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
hold on; plot(squeeze(skilful_CPS_RV(:,floor(window/30))),'k','LineWidth',2); hold off;
xlim([0,70]); ylim([0,1]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_xvar_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
hold on; plot(squeeze(skilful_MRV(:,floor(window/30))),'k','LineWidth',2); hold off;
xlim([0,70]); ylim([0,1]); grid on

% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_xvar_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
hold on; plot(squeeze(skilful_RVM(:,floor(window/30))),'k','LineWidth',2); hold off;
xlim([0,70]); ylim([0,1]); grid on

end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
    hold on; plot([0,70],[skilful_threshold skilful_threshold],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r (',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

% Setting the title
title_str =[];
if strfind(GROUP_NAME,'_nstat') title_str=[title_str, 'NSTAT'];
elseif strfind(GROUP_NAME,'_stat') title_str=[title_str, 'STAT'];
else title_str=[title_str, 'RND'];
end

if strfind(GROUP_NAME,'glb') title_str=[title_str, '_{glb\_']; 
elseif strfind(GROUP_NAME,'ntrop') title_str=[title_str, '_{ntrop\_']; 
else title_str=[title_str, '\_INVALID\_TITLE']; 
end

if strfind(GROUP_NAME,'ts') title_str=[title_str, 'ts}']; 
elseif strfind(GROUP_NAME,'pr') title_str=[title_str, 'pr}']; 
else error('TS/PR if function has failed');
end

suptitle(title_str);


set(gcf, 'PaperPosition', [0 0 20 23]);
legendH = legend('5^t^h Percentile','95^t^h Percentile','Median','location','best','Orientation','horizontal');
set(legendH, 'FontSize',16,'Position',[0.372689213508714 0.0116864099741356 0.54468339307049 0.0419381107491857]);
saveas(gcf,['Plots/calWdwExpVar_vs_NumStns_xvar0.3_',GROUP_NAME,'.jpg']);

% Additional Analysis of rate of skil improvement

% skilful_EPC_RV_30_stat=skilful_EPC_RV;
% skilful_CPS_RV_30_stat=skilful_CPS_RV;
% skilful_MRV_30_stat=skilful_MRV;
% skilful_RVM_30_stat=skilful_RVM;
% 
% skilful_EPC_RV_50_stat=skilful_EPC_RV;
% skilful_CPS_RV_50_stat=skilful_CPS_RV;
% skilful_MRV_50_stat=skilful_MRV;
% skilful_RVM_50_stat=skilful_RVM;
% 
% skilful_EPC_RV_30_nstat=skilful_EPC_RV;
% skilful_CPS_RV_30_nstat=skilful_CPS_RV;
% skilful_MRV_30_nstat=skilful_MRV;
% skilful_RVM_30_nstat=skilful_RVM;
% 
% skilful_EPC_RV_50_nstat=skilful_EPC_RV;
% skilful_CPS_RV_50_nstat=skilful_CPS_RV;
% skilful_MRV_50_nstat=skilful_MRV;
% skilful_RVM_50_nstat=skilful_RVM;

% 0.5 Explained Variance Threshold
clf
letters = 'abcdefghijkl';

s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);

for window=[31 61 91]

    axes(s_Hnd(1+(floor(window/30)-1)*4));
    plot(skilful_EPC_RV_50_stat(:,floor(window/30)) ./ skilful_EPC_RV_50_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(2+(floor(window/30)-1)*4));
    plot(skilful_CPS_RV_50_stat(:,floor(window/30)) ./ skilful_CPS_RV_50_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(3+(floor(window/30)-1)*4));
    plot(skilful_MRV_50_stat(:,floor(window/30)) ./ skilful_MRV_50_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(4+(floor(window/30)-1)*4));
    plot(skilful_RVM_50_stat(:,floor(window/30)) ./ skilful_RVM_50_nstat(:,floor(window/30)),'k','LineWidth',2);
end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    xlim([0,70]); ylim([0,1]); grid on
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r (',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle('Proportion of skilful reconstructions (0.5ExpVar Threshold): STAT_{ntrop\_ts} / NSTAT_{ntrop\_ts}')

% 0.3 Explained Variance Threshold - Divides
clf
letters = 'abcdefghijkl';

s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);

for window=[31 61 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    plot(skilful_EPC_RV_30_stat(:,floor(window/30)) ./ skilful_EPC_RV_30_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(2+(floor(window/30)-1)*4));
    plot(skilful_CPS_RV_30_stat(:,floor(window/30)) ./ skilful_CPS_RV_30_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(3+(floor(window/30)-1)*4));
    plot(skilful_MRV_30_stat(:,floor(window/30)) ./ skilful_MRV_30_nstat(:,floor(window/30)),'k','LineWidth',2);
    axes(s_Hnd(4+(floor(window/30)-1)*4));
    plot(skilful_RVM_30_stat(:,floor(window/30)) ./ skilful_RVM_30_nstat(:,floor(window/30)),'k','LineWidth',2);
end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on','XTick', [0:20:70]); 
    xlim([0,70]); ylim([0,5]); grid on
%     text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(9)); xlabel('Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:1:5]);
%     ylabel(['r (',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle('Proportion of skilful reconstructions (0.3ExpVar Threshold): STAT_{ntrop\_ts} / NSTAT_{ntrop\_ts}')

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
    grid on
    xlim([0 1]); ylim([0 70]); xlabel(' Proportion of non-stationary proxies ','FontSize',16);
    legend([Hnd],'Network Size of 3','Network Size of 5','Network Size of 10','Network Size of 30','Network Size of 50','Network Size of 70');
    % h = findobj(gca,'Type','Patch');
    % set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    ylabel('Percentage of reconstructions','FontSize',16)
%     title('Proportion of nstat stns per reconstruction','FontSize',14);

end

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [0 70], ...
    'XLim'        , [0 1], ...
    'FontSize'    , 16 ...
    );

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [18 14]);
%% Figure 9 - Histogram Comparisons
figure;
window = 31; numstnstocompare=70;
NUM_CAL_WDW = 10; clear CAL_WDW;
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end
s_Hnd = tight_subplot(3,4,[0.02 0.02],[0.1 0.05],[0.1 0.01]);

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

hold off; grid on; xlim([0 1]); ylim([0 0.6]);
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
hold off; grid on; xlim([0 1]); ylim([0 0.6]);
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
hold off; grid on; xlim([0 1]); ylim([0 0.6]);
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
hold off; grid on; xlim([0 1]); ylim([0 0.6]);
% set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick',[0:0.1:0.7],'XTick',[0:0.2:1]); 

i=i+1;
end; clear i;

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    grid on
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on','TickLength'  , [.02 0],'YTick', [0:0.1:0.7],'XTick', [0:0.2:1]); 
    text(0+0.01,0.6-0.05,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(10)); xlabel('Correlation','FontSize',16);
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
set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on','YTick',[0:0.1:0.7]); 
legend([glb_Hnd, ntrop_Hnd],'glb','ntrop','Location','northwest')

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [18 14]);
%% Figure 10
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
subplot(2,1,1);
[c,h]=contourf(lon,lat,squeeze(rc31_eof_ts_fm(1,:,:)),12); plotworld; caxis([-0.03,0.03]); colormap(redblue(12))

set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 2,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );

title('EOF1 (31yr window)','FontSize',16); colorbar;
text(0+1,90-1,'a)','FontSize',22,'FontWeight','bold');
% subplot(3,1,2)
% pcolor(lon,lat,squeeze(rc61_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% title('EOF 2 with window length 61 years')
% subplot(3,1,3)
% pcolor(lon,lat,squeeze(rc91_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% title('EOF 2 with window length 91 years')
% subplot(2,2,4)
subplot(2,2,4)
cmap=hsv(3);
plot(rc31_expvar,'','Color',cmap(1,:),'LineWidth',2); hold on;
plot(rc61_expvar,'','Color',cmap(2,:),'LineWidth',2);
plot(rc91_expvar,'','Color',cmap(3,:),'LineWidth',2); hold off;
ylim([0 50]); xlim([1,10]); legend('31 year window','61 year window','91 year window');
ylabel('Percentage %','FontSize',16); grid on;
xlabel('Number of EOF','FontSize',16)
title('Explained Variance','FontSize',16);
text(0+0.1,1-0.05,'c)','FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );

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
subplot(2,2,3)
a=rc31_PC_ts(1,:)'; b=rc61_PC_ts(1,:)'; c=rc91_PC_ts(1,:)';
plot(15:498-16,a,'','Color',cmap(1,:),'LineWidth',2); hold on;
plot(30:498-31,b,'','Color',cmap(2,:),'LineWidth',2); grid on;
plot(45:498-46,c,'','Color',cmap(3,:),'LineWidth',2); hold off
legend('31 year window','61 year window','91 year window','location','southeast')
xlabel('Year','FontSize',16)
title('Principal Component 1','FontSize',16)
text(0+0.1,1-0.05,'b)','FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );

% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [30 20]);
%% Figure 11
figure;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
for window = [31, 61, 91]

GROUP_NAME = 'pneof_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data-honours/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

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
    text(0+1,1-0.06,[letters(i),')'],'FontSize',20,'FontWeight','bold');
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

% suptitle([strrep(GROUP_NAME,'_','\_')])
suptitle('PNEOF1')
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);

%% Figure 12 - Reconstruction and STD of corr(proxy TS)

% load DataFiles/rvm_vs_mrv.mat;
% first_sig=68; NUM_STNS=10;
% 
% j=1; jj=1; jjj=1; jjjj=1; sig_MRV=[]; sig_RVM=[]; sig_CPS_RV=[]; sig_EPC_RV=[]; 
% for c=1:10
%     for NUM_STNS=10
%         for i=1:NUM_TRIALS
%             recon_MRV = squeeze(MRV_all_grps(c,NUM_STNS,i,16:end-15));
%             recon_RVM = squeeze(RVM_all_grps(c,NUM_STNS,i,16:end-15));
%             recon_CPS_RV = squeeze(CPS_RV_all_grps(c,NUM_STNS,i,16:end-15));
%             recon_EPC_RV = squeeze(EPC_RV_all_grps(c,NUM_STNS,i,16:end-15));
%             std_series = squeeze(std_mov_corr_all_grps(c,NUM_STNS,i,31:end));
%             recon_RVM = (recon_RVM-mean(recon_RVM))./std(recon_RVM);
%             recon_MRV = (recon_MRV-mean(recon_MRV))./std(recon_MRV);
%             recon_CPS_RV = (recon_CPS_RV-mean(recon_CPS_RV))./std(recon_CPS_RV);
%             recon_EPC_RV = (recon_EPC_RV-mean(recon_EPC_RV))./std(recon_EPC_RV);
%             std_series = (std_series-mean(std_series))./std(std_series);
%             [r sig df] = calc_statsig(recon_MRV, std_series,1);
%             [r sig_1 df] = calc_statsig(recon_RVM, std_series,1);
%             [r sig_2 df] = calc_statsig(recon_CPS_RV, std_series,1);
%             [r sig_3 df] = calc_statsig(recon_EPC_RV, std_series,1);
%             
%             if sig>=95
%                 sig_MRV(j,:) = [c, i];
%                 j=j+1;
%             end
%             
%             if sig_1>=95
%                 sig_RVM(jj,:) = [c, i];
%                 jj=jj+1;
%             end
%                 
%             if sig_2>=95 
%                 sig_CPS_RV(jjj,:) = [c, i];
%                 jjj=jjj+1;
%             end
%                 
%             if sig_3>=95
%                 sig_EPC_RV(jjjj,:) = [c, i];
%                 jjjj=jjjj+1;
%             end
%         end
%     end
% end

subplot(3,2,1)
r_meth=squeeze(RVM_all_grps(sig_RVM(first_sig,1),10,sig_RVM(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_RVM(first_sig,1),10,sig_RVM(first_sig,2),31:end));
[axH std_H rvm_H] = plotyy(1:469,temp,1:469,r_meth);
set([std_H rvm_H],'LineWidth',2);
title(['RVM'],'FontSize',16)
text(400,0.05,['r=',num2str(corr(temp,r_meth),'%1.2f')],'FontSize',16)
xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
ylim(axH(1),[0 0.6]); ylim(axH(2),[0 2]);
ylabel(axH(1),'Std','FontSize',16); ylabel(axH(2),'Running Variance','FontSize',16);
set(axH(1),'YTick',[0.1:0.1:0.6],'LineWidth',2,'FontSize',16); set(axH(2),'YTick',[0,1,2],'LineWidth',2,'FontSize',16);
text(0+1,0.6-0.05,['a)'],'FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
subplot(3,2,3)
r_meth=squeeze(CPS_RV_all_grps(sig_CPS_RV(first_sig,1),10,sig_CPS_RV(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_CPS_RV(first_sig,1),10,sig_CPS_RV(first_sig,2),31:end));
[axH std_H rvm_H] = plotyy(1:469,temp,1:469,r_meth);
set([std_H rvm_H],'LineWidth',2);
title(['CPS\_RV'],'FontSize',16)
text(400,0.05,['r=',num2str(corr(temp,r_meth),'%1.2f')],'FontSize',16)
xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
ylim(axH(1),[0,0.6]); ylim(axH(2),[0 2]);
ylabel(axH(1),'Std','FontSize',16); ylabel(axH(2),'Running Variance','FontSize',16);
set(axH(1),'YTick',[0.1:0.1:0.6],'LineWidth',2,'FontSize',16); set(axH(2),'YTick',[0,1,2],'LineWidth',2,'FontSize',16);
text(0+1,0.6-0.05,['c)'],'FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
subplot(3,2,5)
r_meth=squeeze(EPC_RV_all_grps(sig_EPC_RV(first_sig,1),10,sig_EPC_RV(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_EPC_RV(first_sig,1),10,sig_EPC_RV(first_sig,2),31:end));
[axH std_H rvm_H] = plotyy(1:469,temp,1:469,r_meth);
set([std_H rvm_H],'LineWidth',2);
title(['EPC\_RV'],'FontSize',16)
text(400,0.05,['r=',num2str(corr(temp,r_meth),'%1.2f')],'FontSize',16)
xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
ylim(axH(1),[0,0.6]); ylim(axH(2),[0 2]);
ylabel(axH(1),'Std','FontSize',16); ylabel(axH(2),'Running Variance','FontSize',16);
set(axH(1),'YTick',[0.1:0.1:0.6],'LineWidth',2,'FontSize',16); set(axH(2),'YTick',[0,1,2],'LineWidth',2,'FontSize',16);
legend([std_H, rvm_H],'STD of corr(proxy TS)','Reconstruction')
text(0+1,0.6-0.05,['e)'],'FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
% subplot(4,1,4)
% r_meth=squeeze(MRV_all_grps(sig_MRV(first_sig,1),10,sig_MRV(first_sig,2),16:end-15));
% temp=squeeze(std_mov_corr_all_grps(sig_MRV(first_sig,1),10,sig_MRV(first_sig,2),31:end));
% [axH std_H rvm_H] = plotyy(1:469,temp,1:469,r_meth);
% title(['std of corr(proxy temp) vs MRV method, r=',num2str(corr(temp,r_meth))])
% xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
% ylim(axH(1),[0,0.6]); ylim(axH(2),[0 2]);
% ylabel(axH(1),'Std'); ylabel(axH(2),'Running Variance');
% set(axH(1),'YTick',[0:0.1:0.6]); set(axH(2),'YTick',[0,1,2]);
% legend([std_H, rvm_H],'STD of corr(proxy TS)','Reconstruction')


% load DataFiles/rvm_vs_mrv.mat
% i=2; c=1;

subplot(3,2,2)
r_meth=squeeze(RVM_all_grps(sig_RVM(first_sig,1),10,sig_RVM(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_RVM(first_sig,1),10,sig_RVM(first_sig,2),31:end));
plot(temp,r_meth,'.','LineWidth',2); h=lsline; set(h,'LineWidth',2);
title(['RVM'],'FontSize',16)
xlabel('Std of corr(proxy TS)','FontSize',16); ylabel('Reconstruction','FontSize',16);
xlim([0.1 0.35]); ylim([0.25 1.75]);
% xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
% ylim(axH(1),[0 0.6]); ylim(axH(2),[0 2]);
% ylabel(axH(1),'Std'); ylabel(axH(2),'Running Variance');
% set(axH(1),'YTick',[0.1:0.1:0.6]); set(axH(2),'YTick',[0,1,2]);
text(0.1+0.01,1.7-0.03,['b)'],'FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
subplot(3,2,4)
r_meth=squeeze(CPS_RV_all_grps(sig_CPS_RV(first_sig,1),10,sig_CPS_RV(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_CPS_RV(first_sig,1),10,sig_CPS_RV(first_sig,2),31:end));
plot(temp,r_meth,'.','LineWidth',2); h=lsline; set(h,'LineWidth',2);
title(['CPS\_RV'],'FontSize',16)
xlabel('Std of corr(proxy TS)','FontSize',16); ylabel('Reconstruction','FontSize',16);
xlim([0.1 0.35]); ylim([0.25 1.75]);
% xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
% ylim(axH(1),[0,0.6]); ylim(axH(2),[0 2]);
% ylabel(axH(1),'Std'); ylabel(axH(2),'Running Variance');
% set(axH(1),'YTick',[0.1:0.1:0.6]); set(axH(2),'YTick',[0,1,2]);
text(0.1+0.01,1.7-0.03,['d)'],'FontSize',22,'FontWeight','bold');

set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
subplot(3,2,6)
r_meth=squeeze(EPC_RV_all_grps(sig_EPC_RV(first_sig,1),10,sig_EPC_RV(first_sig,2),16:end-15));
temp=squeeze(std_mov_corr_all_grps(sig_EPC_RV(first_sig,1),10,sig_EPC_RV(first_sig,2),31:end));
plot(temp,r_meth,'.','LineWidth',2); h=lsline; set(h,'LineWidth',2);
title(['EPC\_RV'],'FontSize',16)
xlabel('Std of corr(proxy TS)','FontSize',16); ylabel('Reconstruction','FontSize',16);
xlim([0.1 0.35]); ylim([0.25 1.75]);
% xlim(axH(1),[0 469]); xlim(axH(2),[0 469]); grid on;
% ylim(axH(1),[0,0.6]); ylim(axH(2),[0 2]);
% ylabel(axH(1),'Std'); ylabel(axH(2),'Running Variance');
% set(axH(1),'YTick',[0.1:0.1:0.6]); set(axH(2),'YTick',[0,1,2]);
% legend([std_H, rvm_H],'STD of corr(proxy TS)','Reconstruction')
text(0.1+0.01,1.7-0.03,['f)'],'FontSize',22,'FontWeight','bold');


set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'FontSize'    , 16 ...
    );
% % Extra Analysis
% allreg_MRV=nan(size(bad_sig_MRV,1),1);
% for i=1:size(bad_sig_MRV,1)
%     s = squeeze(std_mov_corr_all_grps(bad_sig_MRV(i,1),NUM_STNS,bad_sig_MRV(i,2),(31:end)));
%     r = squeeze(MRV_all_grps(bad_sig_MRV(i,1),NUM_STNS,bad_sig_MRV(i,2),(16:end-15)));
%     [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
%     allreg_MRV(i) = b(1);
% end
% 
% 
% allreg_RVM=nan(size(bad_sig_RVM,1),1);
% for i=1:size(bad_sig_RVM,1)
%     s = squeeze(std_mov_corr_all_grps(bad_sig_RVM(i,1),NUM_STNS,bad_sig_RVM(i,2),(31:end)));
%     r = squeeze(RVM_all_grps(bad_sig_RVM(i,1),NUM_STNS,bad_sig_RVM(i,2),(16:end-15)));
%     [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
%     allreg_RVM(i) = b(1);
% end
% 
% allreg_CPS_RV=nan(size(bad_sig_CPS_RV,1),1);
% for i=1:size(bad_sig_CPS_RV,1)
%     s = squeeze(std_mov_corr_all_grps(bad_sig_CPS_RV(i,1),NUM_STNS,bad_sig_CPS_RV(i,2),(31:end)));
%     r = squeeze(CPS_RV_all_grps(bad_sig_CPS_RV(i,1),NUM_STNS,bad_sig_CPS_RV(i,2),(16:end-15)));
%     [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
%     allreg_CPS_RV(i) = b(1);
% end
% 
% allreg_EPC_RV=nan(size(bad_sig_EPC_RV,1),1);
% for i=1:size(bad_sig_EPC_RV,1)
%     s = squeeze(std_mov_corr_all_grps(bad_sig_EPC_RV(i,1),NUM_STNS,bad_sig_EPC_RV(i,2),(31:end)));
%     r = squeeze(EPC_RV_all_grps(bad_sig_EPC_RV(i,1),NUM_STNS,bad_sig_EPC_RV(i,2),(16:end-15)));
%     [b, b_int, residu, rint, thestats] = regress(r(:),[s(:), ones(size(s(:)))]);
%     allreg_EPC_RV(i) = b(1);
% end
% 
% 
% % Looking at significant bad recons
% a=[]; b=[];
% for c=1:10
%     a = cat(1,a,c*ones(length(bad_recons_MRV{c}(:)),1));
%     b = cat(1,b,bad_recons_MRV{c}(:));
% end
% bad_MRV=[a,b]; 
% bad_sig_MRV=intersect(sig_MRV,[a,b],'rows');
% 
% a=[]; b=[];
% for c=1:10
%     a = cat(1,a,c*ones(length(bad_recons_RVM{c}(:)),1));
%     b = cat(1,b,bad_recons_RVM{c}(:));
% end
% bad_RVM=[a,b]; 
% bad_sig_RVM=intersect(sig_RVM,[a,b],'rows');
% 
% a=[]; b=[];
% for c=1:10
%     a = cat(1,a,c*ones(length(bad_recons_CPS_RV{c}(:)),1));
%     b = cat(1,b,bad_recons_CPS_RV{c}(:));
% end
% bad_CPS_RV=[a,b]; 
% bad_sig_CPS_RV=intersect(sig_CPS_RV,[a,b],'rows');
% 
% a=[]; b=[];
% for c=1:10
%     a = cat(1,a,c*ones(length(bad_recons_EPC_RV{c}(:)),1));
%     b = cat(1,b,bad_recons_EPC_RV{c}(:));
% end
% bad_EPC_RV=[a,b]; 
% bad_sig_EPC_RV=intersect(sig_EPC_RV,[a,b],'rows');
%     

set(gcf, 'PaperSize', [30 20]);
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
folderlist = folderlist([8 11 17 18]+2);% Remove the ones you dont want
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
    contourf(lon,lat,values); plotworld; a=flipud(hot(12)); colormap(a(1:10,:));
    xlim([0 360]); ylim([-90 90]);
    title([strrep(DIR_NAME(58:end),'_','\_'),' - mean: ',num2str(mean(num_prox),'%.0f'),', std: ',num2str(std(num_prox),'%.0f')]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    text(0-5,90-1,[letters(i),')'],'FontSize',20,'FontWeight','bold');       
end

end
% 
% % Proportion of tropical stations
% load('/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/31yrWindow/ntrop_ts/CalWdw:1-31/3stns_1000prox.mat')
% 
% [~,S_bd]= min(abs(lat--10));
% [~,N_bd]= min(abs(lat-10));
% [~,W_bd]= min(abs(lon-100));
% [~,E_bd]= min(abs(lon-300));
% [lat_ind_tmp,lon_ind_tmp]=ind2sub([90 144],indice_pool);
% trop_ind = find(lat_ind_tmp<N_bd & lat_ind_tmp>S_bd & lon_ind_tmp<E_bd & lon_ind_tmp>W_bd);
% 
% length(trop_ind)/length(indice_pool) % Proportion of tropical stations
% % Location Check
% plotworld; hold on; scatter(lon(lon_ind_tmp(trop_ind)),lat(lat_ind_tmp(trop_ind))); hold off;


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
text(0+1,90-1,'b)','FontSize',20,'FontWeight','bold');

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
text(0+1,90-1,'b)','FontSize',20,'FontWeight','bold');

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
text(0+1,90-1,'b)','FontSize',20,'FontWeight','bold');


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
text(0+1,90-1,'b)','FontSize',20,'FontWeight','bold');


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    text(0+1,1-0.06,[letters(i),')'],'FontSize',20,'FontWeight','bold');
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
clf; figure;
letters = 'abcdefghijkl';
numstnstocompare=3:70;
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
for window = [31, 61, 91]

GROUP_NAME = 'pneof_ts'; % Change group name to get other figs
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data-honours/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];

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
    text(0+1,0.5-0.03,[letters(i),')'],'FontSize',20,'FontWeight','bold');
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
saveas(gcf,['../../Dropbox/calWdwRMSE_vs_NumStns_',GROUP_NAME,'.jpg']);


%% Appendix Figure 4


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
temp_var_EPC_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_var_CPS_RV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_var_MRV = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
temp_var_RVM = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'var_EPC_RV','var_CPS_RV','var_MRV','var_RVM')
    temp_var_EPC_RV(:,c,:) = var_EPC_RV;
    temp_var_CPS_RV(:,c,:) = var_CPS_RV;
    temp_var_MRV(:,c,:) = var_MRV;
    temp_var_RVM(:,c,:) = var_RVM;
end

% Plotting EPC

% subplot(3,4,1+(floor(window/30)-1)*4)
axes(s_Hnd(1+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_var_EPC_RV,[.05 .5 .95], 3);
% Range Plotting
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.15]); grid on

% Plotting CPS
% subplot(3,4,2+(floor(window/30)-1)*4)
axes(s_Hnd(2+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_var_CPS_RV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.15]); grid on

% Plotting MRV
% subplot(3,4,3+(floor(window/30)-1)*4)
axes(s_Hnd(3+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_var_MRV,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.15]); grid on


% Plotting RVM
% subplot(3,4,4+(floor(window/30)-1)*4)
axes(s_Hnd(4+(floor(window/30)-1)*4))
corr_RV_qn = quantile(temp_var_RVM,[.05 .5 .95], 3);
corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,1))',squeeze(corr_RV_qn_rng(3:70,1,1))','b','k',[],0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,3))',squeeze(corr_RV_qn_rng(3:70,1,3))','r','k','add',0.5);
jbfill([3:70],squeeze(corr_RV_qn_rng(3:70,2,2))',squeeze(corr_RV_qn_rng(3:70,1,2))','y','k','add',0.5);
xlim([0,70]); ylim([0,0.15]); grid on


end

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',14);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',14);
axes(s_Hnd(3)); title(['MRV'],'FontSize',14);
axes(s_Hnd(4)); title(['RVM'],'FontSize',14);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.05:0.15],'XTick', [0:20:70]);
    text(0+1,0.15-0.01,[letters(i),')'],'FontSize',20,'FontWeight','bold');
end

axes(s_Hnd(9)); xlabel('Proxy Network Size');
for window = [31, 61, 91]
    axes(s_Hnd(1+(floor(window/30)-1)*4));
    set(gca,'YTickLabel',[0:0.05:0.15]);
    ylabel(['sigma^2 (',num2str(window),'yrs)'])
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

% suptitle([strrep(GROUP_NAME,'_','\_')])
suptitle('RND_{glb}');
set(gcf, 'PaperPosition', [0 0 19 23]);
legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','best','Orientation','horizontal');
set(legendH, 'FontSize',10);
