% This script plots the new figures in the CP Paper from the results
% section onwards.

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


%% Figure Alpha
figure;

letters = 'abcdefghijkl';
s_Hnd = tight_subplot(2,2,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

load('DataFiles/500yrCalWdw_meth_stats.mat')
qEPC = quantile(all_stn_corr_EPC_RV,[.05 .5 .95], 2);
qMRV = quantile(all_stn_corr_MRV,[.05 .5 .95], 2);
qCPS = quantile(all_stn_corr_CPS_RV,[.05 .5 .95], 2);
qRVM = quantile(all_stn_corr_RVM,[.05 .5 .95], 2);

axes(s_Hnd(1)); hold on;
Hnd(1,1) = plot(squeeze(qEPC(:,1)),'--k','LineWidth',2);
Hnd(1,3) = plot(squeeze(qEPC(:,3)),'-.k','LineWidth',2);
Hnd(1,2) = plot(squeeze(qEPC(:,2)),'k','LineWidth',3);

axes(s_Hnd(2)); hold on;
Hnd(2,1) = plot(squeeze(qCPS(:,1)),'--k','LineWidth',2);
Hnd(2,3) = plot(squeeze(qCPS(:,3)),'-.k','LineWidth',2);
Hnd(2,2) = plot(squeeze(qCPS(:,2)),'k','LineWidth',3);

axes(s_Hnd(3)); hold on;
Hnd(3,1) = plot(squeeze(qMRV(:,1)),'--k','LineWidth',2);
Hnd(3,3) = plot(squeeze(qMRV(:,3)),'-.k','LineWidth',2);
Hnd(3,2) = plot(squeeze(qMRV(:,2)),'k','LineWidth',3);

axes(s_Hnd(4)); hold on;
Hnd(4,1) = plot(squeeze(qRVM(:,1)),'--k','LineWidth',2);
Hnd(4,3) = plot(squeeze(qRVM(:,3)),'-.k','LineWidth',2);
Hnd(4,2) = plot(squeeze(qRVM(:,2)),'k','LineWidth',3);

window = 31;

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
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

% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w31_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);


window = 61; clear CAL_WDW;
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_xvar_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_xvar_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_xvar_MRV(:,c,:) = all_stn_corr_MRV;
    temp_xvar_RVM(:,c,:) = all_stn_corr_RVM;
end
    
% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w61_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
    

window = 91; clear CAL_WDW;
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_xvar_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_xvar_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_xvar_MRV(:,c,:) = all_stn_corr_MRV;
    temp_xvar_RVM(:,c,:) = all_stn_corr_RVM;
end
    
% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w91_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);


% Formatting
axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
legend(s_Hnd(1),[Hnd(1,2) Hnd(1,1) Hnd(1,3)],'Median (499yr)', ...
       '5^t^h %ile (499yr) ','95^t^h %ile (499yr)','location','southeast'                );
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
legend(s_Hnd(2),[w91_Hnd(1,2) w91_Hnd(1,1) w91_Hnd(1,3)],'Median (91yr)', ...
       '5^t^h %ile (91yr) ','95^t^h %ile (91yr)','location','southeast'                );
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
legend(s_Hnd(3),[w61_Hnd(1,2) w61_Hnd(1,1) w61_Hnd(1,3)],'Median (61yr)', ...
       '5^t^h %ile (61yr) ','95^t^h %ile (61yr)','location','southeast'                );
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);
legend(s_Hnd(4),[w31_Hnd(1,2) w31_Hnd(1,1) w31_Hnd(1,3)],'Median (31yr)', ...
       '5^t^h %ile (31yr) ','95^t^h %ile (31yr)','location','southeast'                );

for i=1:4
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    xlim([0 70]); ylim([0.2 1]); grid on
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(3)); xlabel('Network Size');
for i=1:2
    axes(s_Hnd(i*2-1));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r'])
end

for i=3:4
    axes(s_Hnd(i));
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



%% Figure Bravo
clf
letters = 'abcdefghijkl';

% Skilful Threshold
skilful_threshold = sqrt(0.5); % In correlation

s_Hnd = tight_subplot(3,4,[0.01 0.005],[0.10 0.01],[0.1 0.01]);
skilful_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
skilful_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
skilful_MRV = nan(max(numstnstocompare),length([31 61 91]));
skilful_RVM = nan(max(numstnstocompare),length([31 61 91]));
w_Hnd = nan(3,3);
GROUP_NAME = {'glb_ts','ntrop_ts','ntrop_ts_nstat'}; % Change group name to get other figs

for group_name_counter = 1:3 
    
for window = [31, 61, 91]

    % DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
    % DIR_NAME = ['/home/nfs/z3372730/Documents/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',GROUP_NAME{group_name_counter}];
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
    axes(s_Hnd(1+(group_name_counter-1)*4)); hold on;
    corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
    if window == 31
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','b','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[0.6 0.6 1]);
    elseif window == 61
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','y','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
% %         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 1]);
    elseif window == 91
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 0.6]);
%         Leg_H = legend(s_Hnd(1),[w_Hnd(1,2) w_Hnd(1,1) w_Hnd(1,3)],'Median (31yr)', ...
%         '5^t^h %ile (31yr) ','95^t^h %ile (31yr)','location','southeast');
%         set(Leg_H, 'FontSize',14);
    end
        
    % Plotting CPS
    axes(s_Hnd(2+(group_name_counter-1)*4)); hold on;
    corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
    if window == 31
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','b','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[0.6 0.6 1]);
    elseif window == 61
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','y','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 1]);
    elseif window == 91
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 0.6]);
%         Leg_H = legend(s_Hnd(2),[w_Hnd(2,2) w_Hnd(2,1) w_Hnd(2,3)],'Median (61yr)', ...
%         '5^t^h %ile (61yr) ','95^t^h %ile (61yr)','location','southeast');
%         set(Leg_H, 'FontSize',14);
    end
    
    % Plotting MRV
    axes(s_Hnd(3+(group_name_counter-1)*4)); hold on;
    corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
    if window == 31
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','b','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[0.6 0.6 1]);
    elseif window == 61
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','y','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 1]);
    elseif window == 91
        jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 0.6]);
%         Leg_H = legend(s_Hnd(3),[w_Hnd(3,2) w_Hnd(3,1) w_Hnd(3,3)],'Median (91yr)', ...
%         '5^t^h %ile (91yr) ','95^t^h %ile (91yr)','location','southeast');
%         set(Leg_H, 'FontSize',14);
    end
    
    % Plotting RVM
    axes(s_Hnd(4+(group_name_counter-1)*4)); hold on;
    corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
    if window == 31
        f_Hnd(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','b','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[0.6 0.6 1]);
    elseif window == 61
        f_Hnd(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','y','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 1]);
    elseif window == 91
        f_Hnd(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
%         w_Hnd(floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
%         w_Hnd(floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
%         w_Hnd(floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);
%         set(w_Hnd(floor(window/30),:),'Color',[1 0.6 0.6]);
    end
    
    if window == 31
        axes(s_Hnd(1+(group_name_counter-1)*4));
        h31(1)=plot(skilful_EPC_RV(:,floor(window/30)),'b','LineWidth',2);
        axes(s_Hnd(2+(group_name_counter-1)*4));
        h31(2)=plot(skilful_CPS_RV(:,floor(window/30)),'b','LineWidth',2);
        axes(s_Hnd(3+(group_name_counter-1)*4));
        h31(3)=plot(skilful_MRV(:,floor(window/30)),'b','LineWidth',2);
        axes(s_Hnd(4+(group_name_counter-1)*4));
        h31(4)=plot(skilful_RVM(:,floor(window/30)),'b','LineWidth',2);

    elseif window == 61
        axes(s_Hnd(1+(group_name_counter-1)*4));
        h61(1)=plot(skilful_EPC_RV(:,floor(window/30)),'y','LineWidth',2);
        axes(s_Hnd(2+(group_name_counter-1)*4));
        h61(2)=plot(skilful_CPS_RV(:,floor(window/30)),'y','LineWidth',2);
        axes(s_Hnd(3+(group_name_counter-1)*4));
        h61(3)=plot(skilful_MRV(:,floor(window/30)),'y','LineWidth',2);
        axes(s_Hnd(4+(group_name_counter-1)*4));
        h61(4)=plot(skilful_RVM(:,floor(window/30)),'y','LineWidth',2);
    elseif window == 91
        axes(s_Hnd(1+(group_name_counter-1)*4));
        h91(1)=plot(skilful_EPC_RV(:,floor(window/30)),'r','LineWidth',2);
        axes(s_Hnd(2+(group_name_counter-1)*4));
        h91(2)=plot(skilful_CPS_RV(:,floor(window/30)),'r','LineWidth',2);
        axes(s_Hnd(3+(group_name_counter-1)*4));
        h91(3)=plot(skilful_MRV(:,floor(window/30)),'r','LineWidth',2);
        axes(s_Hnd(4+(group_name_counter-1)*4));
        h91(4)=plot(skilful_RVM(:,floor(window/30)),'r','LineWidth',2);

    end
       
end

for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
    uistack(h31(i),'top');
    uistack(h61(i),'top');
    uistack(h91(i),'top');
end

end

Leg_H = legend(s_Hnd(7),[h31(1) h61(1) h91(1)],'Proportion (31yr)', ...
        'Proportion (61yr) ','Proportion (91yr)','location','southeast');
set(Leg_H, 'FontSize',14);
f_leg_H = legend(s_Hnd(8),[f_Hnd(1) f_Hnd(2) f_Hnd(3)],'%ile range (31yr)', ... 
                 '%ile range (61yr)','%ile range (91yr)','location','southeast');
set(f_leg_H, 'FontSize',14);

axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);

for i=1:12
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    xlim([0 70]); ylim([0 1]); grid on;
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
    hold on; st_H = plot([0,70],[skilful_threshold skilful_threshold],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
    uistack(st_H,'down');
end

st_LH = legend(s_Hnd(6),st_H, 'Skill Threshold');
set(st_LH, 'FontSize', 14);

axes(s_Hnd(9)); xlabel('Network Size');
for group_name_counter=1:3
    axes(s_Hnd(1+(group_name_counter-1)*4));
    set(gca,'YTickLabel',[0:0.2:1]);
    
    % Setting the title
    title_str =[];
    if strfind(GROUP_NAME{group_name_counter},'_nstat') title_str=[title_str, 'NSTAT'];
    elseif strfind(GROUP_NAME{group_name_counter},'_stat') title_str=[title_str, 'STAT'];
    else title_str=[title_str, 'RND'];
    end

    if strfind(GROUP_NAME{group_name_counter},'glb') title_str=[title_str, '_{glb\_']; 
    elseif strfind(GROUP_NAME{group_name_counter},'ntrop') title_str=[title_str, '_{ntrop\_']; 
    else title_str=[title_str, '\_INVALID\_TITLE']; 
    end

    if strfind(GROUP_NAME{group_name_counter},'ts') title_str=[title_str, 'ts}']; 
    elseif strfind(GROUP_NAME{group_name_counter},'pr') title_str=[title_str, 'pr}']; 
    else error('TS/PR if function has failed');
    end
    
    ylabel([title_str])

end


for i=1:4
    axes(s_Hnd(i+8));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle('Window Length and Reconstruction Skill');


%% Figure Charlie

clf
letters = 'abcdefghijkl';
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);

% Skilful Threshold
skilful_threshold = sqrt(0.5); % In correlation

skilful_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
skilful_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
skilful_MRV = nan(max(numstnstocompare),length([31 61 91]));
skilful_RVM = nan(max(numstnstocompare),length([31 61 91]));
w_Hnd_stat = nan(4,3,3); w_Hnd_nstat = nan(4,3,3); prop_Hnd = nan(4,3,3);

for window = [31, 61, 91]

GROUP_NAME = 'ntrop_ts_stat'; % Change group name to get other figs
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
axes(s_Hnd(1+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','g','k',[],0.3);
% w_Hnd_stat(1,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(1,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(1,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'g','Color',[0.5 1 0.5],'LineWidth',2);

% Plotting CPS
axes(s_Hnd(2+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','g','k',[],0.3);
% w_Hnd_stat(2,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(2,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(2,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'g','Color',[0.5 1 0.5],'LineWidth',2);

% Plotting MRV
axes(s_Hnd(3+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','g','k',[],0.3);
% w_Hnd_stat(3,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(3,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(3,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'g','Color',[0.5 1 0.5],'LineWidth',2);

% Plotting RVM
axes(s_Hnd(4+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','g','k',[],0.3);
% w_Hnd_stat(4,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(4,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.g','Color',[0.5 1 0.5],'LineWidth',1);
% w_Hnd_stat(4,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'g','Color',[0.5 1 0.5],'LineWidth',2);

end

% Additional Analysis of rate of skil improvement

skilful_EPC_RV_50_stat=skilful_EPC_RV;
skilful_CPS_RV_50_stat=skilful_CPS_RV;
skilful_MRV_50_stat=skilful_MRV;
skilful_RVM_50_stat=skilful_RVM;


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
axes(s_Hnd(1+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
% w_Hnd_nstat(1,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(1,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(1,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);

% Plotting CPS
axes(s_Hnd(2+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
% w_Hnd_nstat(2,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(2,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(2,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);

% Plotting MRV
axes(s_Hnd(3+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
% w_Hnd_nstat(3,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(3,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(3,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);

% Plotting RVM
axes(s_Hnd(4+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))','r','k',[],0.3);
% w_Hnd_nstat(4,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(4,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(4,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);


end

skilful_EPC_RV_50_nstat=skilful_EPC_RV;
skilful_CPS_RV_50_nstat=skilful_CPS_RV;
skilful_MRV_50_nstat=skilful_MRV;
skilful_RVM_50_nstat=skilful_RVM;


% ***********************0.5 Explained Variance Threshold**********************

for window=[31 61 91]

    axes(s_Hnd(1+(floor(window/30)-1)*4));
    prop_Hnd(1,floor(window/30),1) = plot(skilful_EPC_RV_50_stat(:,floor(window/30)),'Color',[0 0.8 0],'LineWidth',3); 
    prop_Hnd(1,floor(window/30),2) = plot(skilful_EPC_RV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(1,floor(window/30),3) = plot(skilful_EPC_RV_50_stat(:,floor(window/30)) - skilful_EPC_RV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(2+(floor(window/30)-1)*4));
    prop_Hnd(2,floor(window/30),1) = plot(skilful_CPS_RV_50_stat(:,floor(window/30)),'Color',[0 0.8 0],'LineWidth',3);  
    prop_Hnd(2,floor(window/30),2) = plot(skilful_CPS_RV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(2,floor(window/30),3) = plot(skilful_CPS_RV_50_stat(:,floor(window/30)) - skilful_CPS_RV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(3+(floor(window/30)-1)*4));
    prop_Hnd(3,floor(window/30),1) = plot(skilful_MRV_50_stat(:,floor(window/30)),'Color',[0 0.8 0],'LineWidth',3); 
    prop_Hnd(3,floor(window/30),2) = plot(skilful_MRV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(3,floor(window/30),3) = plot(skilful_MRV_50_stat(:,floor(window/30)) - skilful_MRV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(4+(floor(window/30)-1)*4));
    prop_Hnd(4,floor(window/30),1) = plot(skilful_RVM_50_stat(:,floor(window/30)),'Color',[0 0.8 0],'LineWidth',3);  
    prop_Hnd(4,floor(window/30),2) = plot(skilful_RVM_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3);  
    prop_Hnd(4,floor(window/30),3) = plot(skilful_RVM_50_stat(:,floor(window/30)) - skilful_RVM_50_nstat(:,floor(window/30)),'k','LineWidth',3);
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
    hold on; st_H = plot([0,70],[skilful_threshold skilful_threshold],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
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

% legend(s_Hnd(3),[f_Hnd_stat(3,1,2) f_Hnd_stat(3,1,1) f_Hnd_stat(3,1,3)],'stat Median','5^t^h %ile','95^t^h %ile','location','southeast');
% legend(s_Hnd(7),[f_Hnd_nstat(3,1,2) f_Hnd_nstat(3,1,1) f_Hnd_nstat(3,1,3)],'nstat Median','5^t^h %ile','95^t^h %ile','location','southeast');
% legend(s_Hnd(11),[prop_Hnd(3,1,1) prop_Hnd(3,1,2) prop_Hnd(3,1,3)],'Skilful stat','Skilful nstat','Difference');
% st_LH = legend(s_Hnd(6),st_H, 'Skill Threshold');
big_H = legend(s_Hnd(11),[f_Hnd_stat(1),f_Hnd_nstat(1),prop_Hnd(3,1,1) prop_Hnd(3,1,2) prop_Hnd(3,1,3),st_H],...
       'Stationary','Non-stationary','Skilful stat','Skilful nstat','Difference','Skill Threshold','orientation','horizontal');
set(big_H, 'FontSize', 14);

suptitle('STAT_{ntrop\_ts} - NSTAT_{ntrop\_ts}')


%% Figure Delta

figure;

letters = 'abcdefghijkl';
s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

window = 31;

GROUP_NAME = 'pneof_ts'; % Change group name to get other figs
% DIR_NAME = ['../Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
% DIR_NAME = ['/home/nfs/z3372730/Documents/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data-honours/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
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

% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w31_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w31_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--b','LineWidth',1);
w31_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.b','LineWidth',1);
w31_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'b','LineWidth',2);


window = 61; clear CAL_WDW;
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data-honours/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_xvar_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_xvar_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_xvar_MRV(:,c,:) = all_stn_corr_MRV;
    temp_xvar_RVM(:,c,:) = all_stn_corr_RVM;
end
    
% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w61_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w61_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--m','LineWidth',1);
w61_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.m','LineWidth',1);
w61_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'m','LineWidth',2);
    

window = 91; clear CAL_WDW;
DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data-honours/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
end

for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
     'all_stn_corr_EPC_RV','all_stn_corr_CPS_RV','all_stn_corr_MRV','all_stn_corr_RVM')
    temp_xvar_EPC_RV(:,c,:) = all_stn_corr_EPC_RV;
    temp_xvar_CPS_RV(:,c,:) = all_stn_corr_CPS_RV;
    temp_xvar_MRV(:,c,:) = all_stn_corr_MRV;
    temp_xvar_RVM(:,c,:) = all_stn_corr_RVM;
end
    
% Plotting EPC
axes(s_Hnd(1)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
w91_Hnd(1,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(1,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(1,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting CPS
axes(s_Hnd(2)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting MRV
axes(s_Hnd(3)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);

% Plotting RVM
axes(s_Hnd(4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
w91_Hnd(2,1) = plot(squeeze(corr_RV_qn(:,1)),'--r','LineWidth',1);
w91_Hnd(2,3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','LineWidth',1);
w91_Hnd(2,2) = plot(squeeze(corr_RV_qn(:,2)),'r','LineWidth',2);


% Formatting
axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
leg_h = legend(s_Hnd(2),[w91_Hnd(1,2) w91_Hnd(1,3) w91_Hnd(1,1)],'Median (91yr)', ...
       '5^t^h %ile (91yr) ','95^t^h %ile (91yr)','location','southeast'                );
set(leg_h,'FontSize',14,'orientation','horizontal');
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
leg_h = legend(s_Hnd(3),[w61_Hnd(1,2) w61_Hnd(1,3) w61_Hnd(1,1)],'Median (61yr)', ...
       '5^t^h %ile (61yr) ','95^t^h %ile (61yr)','location','southeast'                );
set(leg_h,'FontSize',14,'orientation','horizontal');
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);
leg_h = legend(s_Hnd(4),[w31_Hnd(1,2) w31_Hnd(1,3) w31_Hnd(1,1)],'Median (31yr)', ...
       '5^t^h %ile (31yr) ','95^t^h %ile (31yr)','location','southeast'                );
set(leg_h,'FontSize',14,'orientation','horizontal');

for i=1:4
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    xlim([0 70]); ylim([0 1]); grid on
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(1)); xlabel('Network Size');
for i=1:1
    axes(s_Hnd(i*2-1));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['r'])
end

for i=1:4
    axes(s_Hnd(i));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle('PNEOF1');

%*********************Just manually Stitch together************************
% EOF Maps
% load('DataFiles/all_rcor_ts_eofs.mat')
% subplot(2,1,1);
% [c,h]=contourf(lon,lat,squeeze(rc31_eof_ts_fm(1,:,:)),12); plotworld; caxis([-0.03,0.03]); colormap(redblue(12))
% 
% set(gca, ...
%     'TickDir','out', ...
%     'Box', 'on',    ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'off'      , ...
%     'YTick'       , -90:30:90, ...
%     'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
%     'LineWidth'   , 2,  ...    
%     'YLim'        , [-90 90], ...
%     'XLim'        , [0 360],  ...
%     'FontSize'    , 16   ...
%     );
% 
% title('EOF1 (31yr window)','FontSize',16); colorbar;
% text(0+1,90-1,'a)','FontSize',22,'FontWeight','bold');
% % subplot(3,1,2)
% % pcolor(lon,lat,squeeze(rc61_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% % title('EOF 2 with window length 61 years')
% % subplot(3,1,3)
% % pcolor(lon,lat,squeeze(rc91_eof_ts_fm(2,:,:))); plotworld; colormap(b2r(-0.03,0.03))
% % title('EOF 2 with window length 91 years')
% % subplot(2,2,4)
% subplot(2,2,4)
% cmap=hsv(3);
% plot(rc31_expvar,'','Color',cmap(1,:),'LineWidth',2); hold on;
% plot(rc61_expvar,'','Color',cmap(2,:),'LineWidth',2);
% plot(rc91_expvar,'','Color',cmap(3,:),'LineWidth',2); hold off;
% ylim([0 50]); xlim([1,10]); legend('31 year window','61 year window','91 year window');
% ylabel('Percentage %','FontSize',16); grid on;
% xlabel('Number of EOF','FontSize',16)
% title('Explained Variance','FontSize',16);
% text(0+0.1,1-0.05,'c)','FontSize',22,'FontWeight','bold');
% 
% set(gca, ...
%     'TickDir','in', ...
%     'Box', 'on',    ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'off'      , ...
%     'YMinorTick'  , 'off'      , ...
%     'LineWidth'   , 2,  ...    
%     'FontSize'    , 16 ...
%     );
% 
% % Correlation between PC timeseries
% % a=rc61_PC_ts(2,:)'; b=rc91_PC_ts(2,:)';
% % a=rc61_eof_ts_fm(2,:)'; b=rc91_eof_ts_fm(2,:)';
% % corr(a,b)
% % cor = 0;
% % for i=0:30
% %     cor(i+1) = corr(a((i+1):400+i),b(1:400))
% % end
% % plot(cor); grid minor
% 
% % Plotting PC Time series
% subplot(2,2,3)
% a=rc31_PC_ts(1,:)'; b=rc61_PC_ts(1,:)'; c=rc91_PC_ts(1,:)';
% plot(15:498-16,a,'','Color',cmap(1,:),'LineWidth',2); hold on;
% plot(30:498-31,b,'','Color',cmap(2,:),'LineWidth',2); grid on;
% plot(45:498-46,c,'','Color',cmap(3,:),'LineWidth',2); hold off
% legend('31 year window','61 year window','91 year window','location','southeast')
% xlabel('Year','FontSize',16)
% title('Principal Component 1','FontSize',16)
% text(0+0.1,1-0.05,'b)','FontSize',22,'FontWeight','bold');
% 
% set(gca, ...
%     'TickDir','in', ...
%     'Box', 'on',    ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'off'      , ...
%     'YMinorTick'  , 'off'      , ...
%     'LineWidth'   , 2,  ...    
%     'FontSize'    , 16 ...
%     );
% 
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
% set(gcf, 'PaperSize', [30 20]);

%% Figure Echo

% Existing Figure 8

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

%% Figure Foxtrot

clf
letters = 'abcdefghijkl';
s_Hnd = tight_subplot(3,4,[0.01 0.01],[0.10 0.01],[0.1 0.01]);
glb_colour = [0.5 0.5 1];
ntrop_colour = [1 0.5 0.5];
shading_opacity = 0.5;

% Skilful Threshold
skilful_threshold = sqrt(0.5); % In correlation

skilful_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
skilful_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
skilful_MRV = nan(max(numstnstocompare),length([31 61 91]));
skilful_RVM = nan(max(numstnstocompare),length([31 61 91]));
w_Hnd_stat = nan(4,3,3); w_Hnd_nstat = nan(4,3,3); prop_Hnd = nan(4,3,3);

for window = [31, 61, 91]

GROUP_NAME = 'glb_ts'; % Change group name to get other figs
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
axes(s_Hnd(1+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',glb_colour,'k',[],shading_opacity);

% Plotting CPS
axes(s_Hnd(2+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',glb_colour,'k',[],shading_opacity);

% Plotting MRV
axes(s_Hnd(3+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',glb_colour,'k',[],shading_opacity);

% Plotting RVM
axes(s_Hnd(4+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
f_Hnd_stat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',glb_colour,'k',[],shading_opacity);

end

% Additional Analysis of rate of skil improvement

skilful_EPC_RV_50_stat=skilful_EPC_RV;
skilful_CPS_RV_50_stat=skilful_CPS_RV;
skilful_MRV_50_stat=skilful_MRV;
skilful_RVM_50_stat=skilful_RVM;


skilful_EPC_RV = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 499yr recons
skilful_CPS_RV = nan(max(numstnstocompare),length([31 61 91]));
skilful_MRV = nan(max(numstnstocompare),length([31 61 91]));
skilful_RVM = nan(max(numstnstocompare),length([31 61 91]));

for window = [31, 61, 91]

GROUP_NAME = 'ntrop_ts'; % Change group name to get other figs
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

% % Plotting EPC
% axes(s_Hnd(1+(floor(window/30)-1)*4)); hold on;
% corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
% w_Hnd_nstat(1,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(1,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(1,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);
% 
% % Plotting CPS
% axes(s_Hnd(2+(floor(window/30)-1)*4)); hold on;
% corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
% w_Hnd_nstat(2,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(2,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(2,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);
% 
% % Plotting MRV
% axes(s_Hnd(3+(floor(window/30)-1)*4)); hold on;
% corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
% w_Hnd_nstat(3,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(3,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(3,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);
% 
% % Plotting RVM
% axes(s_Hnd(4+(floor(window/30)-1)*4)); hold on;
% corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
% w_Hnd_nstat(4,floor(window/30),1) = plot(squeeze(corr_RV_qn(:,1)),'--r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(4,floor(window/30),3) = plot(squeeze(corr_RV_qn(:,3)),'-.r','Color',[1 0.5 0.5],'LineWidth',1);
% w_Hnd_nstat(4,floor(window/30),2) = plot(squeeze(corr_RV_qn(:,2)),'r','Color',[1 0.5 0.5],'LineWidth',2);


% Plotting EPC
axes(s_Hnd(1+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_EPC_RV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',ntrop_colour,'k',[],shading_opacity);

% Plotting CPS
axes(s_Hnd(2+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_CPS_RV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',ntrop_colour,'k',[],shading_opacity);

% Plotting MRV
axes(s_Hnd(3+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_MRV(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',ntrop_colour,'k',[],shading_opacity);

% Plotting RVM
axes(s_Hnd(4+(floor(window/30)-1)*4)); hold on;
corr_RV_qn = quantile(temp_xvar_RVM(:,:),[.05 .5 .95], 2);
f_Hnd_nstat(floor(window/30)) = jbfill([3:70],squeeze(corr_RV_qn(3:70,1))',squeeze(corr_RV_qn(3:70,3))',ntrop_colour,'k',[],shading_opacity);

end

skilful_EPC_RV_50_nstat=skilful_EPC_RV;
skilful_CPS_RV_50_nstat=skilful_CPS_RV;
skilful_MRV_50_nstat=skilful_MRV;
skilful_RVM_50_nstat=skilful_RVM;


% ***********************0.5 Explained Variance Threshold**********************

for window=[31 61 91]

    axes(s_Hnd(1+(floor(window/30)-1)*4));
    prop_Hnd(1,floor(window/30),1) = plot(skilful_EPC_RV_50_stat(:,floor(window/30)),'Color',[0 0 0.8],'LineWidth',3); 
    prop_Hnd(1,floor(window/30),2) = plot(skilful_EPC_RV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(1,floor(window/30),3) = plot(skilful_EPC_RV_50_stat(:,floor(window/30)) - skilful_EPC_RV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(2+(floor(window/30)-1)*4));
    prop_Hnd(2,floor(window/30),1) = plot(skilful_CPS_RV_50_stat(:,floor(window/30)),'Color',[0 0 0.8],'LineWidth',3);  
    prop_Hnd(2,floor(window/30),2) = plot(skilful_CPS_RV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(2,floor(window/30),3) = plot(skilful_CPS_RV_50_stat(:,floor(window/30)) - skilful_CPS_RV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(3+(floor(window/30)-1)*4));
    prop_Hnd(3,floor(window/30),1) = plot(skilful_MRV_50_stat(:,floor(window/30)),'Color',[0 0 0.8],'LineWidth',3); 
    prop_Hnd(3,floor(window/30),2) = plot(skilful_MRV_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3); 
    prop_Hnd(3,floor(window/30),3) = plot(skilful_MRV_50_stat(:,floor(window/30)) - skilful_MRV_50_nstat(:,floor(window/30)),'k','LineWidth',3);
    axes(s_Hnd(4+(floor(window/30)-1)*4));
    prop_Hnd(4,floor(window/30),1) = plot(skilful_RVM_50_stat(:,floor(window/30)),'Color',[0 0 0.8],'LineWidth',3);  
    prop_Hnd(4,floor(window/30),2) = plot(skilful_RVM_50_nstat(:,floor(window/30)),'Color',[0.8 0 0],'LineWidth',3);  
    prop_Hnd(4,floor(window/30),3) = plot(skilful_RVM_50_stat(:,floor(window/30)) - skilful_RVM_50_nstat(:,floor(window/30)),'k','LineWidth',3);
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
    hold on; st_H = plot([0,70],[skilful_threshold skilful_threshold],'Color',[0.5 0.5 0.5],'LineWidth',1.5); hold off;
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

big_H = legend(s_Hnd(11),[f_Hnd_stat(1),f_Hnd_nstat(1),prop_Hnd(3,1,1) prop_Hnd(3,1,2) prop_Hnd(3,1,3),st_H],...
       'global','non-tropical','Skilful glb','Skilful ntrop','Difference','Skill Threshold','orientation','horizontal');
set(big_H, 'FontSize', 14);

suptitle('RND_{glb\_ts} - RND_{ntrop\_ts}')
