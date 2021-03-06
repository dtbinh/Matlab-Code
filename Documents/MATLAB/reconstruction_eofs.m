% This script will compute the first few EOFs of the running correlations of the model data and show the
% variance explained
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

%% Model EOFs

NUM_OF_EOFS = 5;
ats_fm = reshape(ats,size(ats,1),size(ats,2)*size(ats,3));
[eof_ats,PC_ats,expvar_ats] = caleof(ats_fm, NUM_OF_EOFS, 2);
eof_ats_fm = reshape(eof_ats,NUM_OF_EOFS,size(ats,2),size(ats,3));
display('EOF_ats done');
%% Finding EOFs of reconstructions

% More Setup
GROUP_NAME = 'glb_ts_nstat_sigpcd';
NUM_OF_EOFS = 5; NUM_TRIALS = 1000; NUM_YRS = 499; numstnstocompare = [3:70]; NUM_CAL_WDW = 10;

all_eof_EPC = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_TRIALS,'single');
all_PC_EPC = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_YRS,'single');
all_expvar_EPC = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS,'single');
all_eof_CPS = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_TRIALS,'single');
all_PC_CPS = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_YRS,'single');
all_expvar_CPS = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS,'single');
corr_all_PC_EPC = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_OF_EOFS);
corr_all_PC_CPS = nan(3,numstnstocompare(end), NUM_CAL_WDW, NUM_OF_EOFS, NUM_OF_EOFS);
tic;
for window = [31, 61, 91]
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
    CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
    end
    DIR_NAME = ['/srv/ccrc/data34/z3372730/Katana_Data/Data/Pseudoproxies/',num2str(window),'yrWindow/',num2str(GROUP_NAME)];
for c=1:size(CAL_WDW,1)
    load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],'all_stn_CPS','all_stn_EPC');
    
    for n=numstnstocompare
            clear eof_EPC PC_EPC expvar_EPC eof_CPS PC_CPS expvar_CPS
            [eof_EPC,PC_EPC,expvar_EPC] = caleof(squeeze(double(all_stn_EPC(n,:,:)))', NUM_OF_EOFS, 2);
            [eof_CPS,PC_CPS,expvar_CPS] = caleof(squeeze(double(all_stn_CPS(n,:,:)))', NUM_OF_EOFS, 2);
            corr_all_PC_EPC(floor(window/30),n,c,:,:) = corr(PC_EPC',PC_ats');
            corr_all_PC_CPS(floor(window/30),n,c,:,:) = corr(PC_CPS',PC_ats'); % first arg becomes row, second is column in corr matrix
            
            all_eof_EPC(floor(window/30),n,c,:,:) = eof_EPC;
            all_PC_EPC(floor(window/30),n,c,:,:) = PC_EPC;
            all_expvar_EPC(floor(window/30),n,c,:,:) = expvar_EPC;
            all_eof_CPS(floor(window/30),n,c,:,:) = eof_CPS;
            all_PC_CPS(floor(window/30),n,c,:,:) = PC_CPS;
            all_expvar_CPS(floor(window/30),n,c,:,:) = expvar_CPS;

            toc; % Took 1200s in total for all
    end
    
    
end
end

save(['DataFiles/reconEOFs_',GROUP_NAME,'.mat'],'all_eof_EPC','all_PC_EPC','all_expvar_EPC', ...
     'all_eof_CPS','all_PC_CPS','all_expvar_CPS', 'corr_all_PC_EPC','corr_all_PC_CPS');

%% Plotting yay

GROUP_NAME = 'glb_ts';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rcorwdw=1; groups=3:70; calwdw=1; recEOFs=1:5; modEOFs=1:5; 

% EOF vs EOF corr plot
axis equal;
contourf(squeeze(mean(mean(mean(abs(corr_all_PC_EPC(rcorwdw,groups,calwdw,recEOFs,modEOFs)),1),2),3)),10);
grid on; colorbar; caxis([0,1]); colormap(flipud(gray(10))); % The rows of contourf plot are same rows in corr  matrix
xlabel('EOF of Model temperature'); ylabel('EOF of nino3.4 Reconstruction')
set(gca,  'XTick', [1:5], 'YTick', [1:5]); 
title(['Mean abs corr for ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs), groups=',...
      num2str(groups(1)),':',num2str(groups(end)),', calwdw=',...
      num2str(calwdw(1)),':',num2str(calwdw(end))                                    ])
  
set(gcf, 'PaperUnits', 'centimeters'); % May already be default
    set(gcf, 'PaperPosition', [0 0 19 17]);
saveas(gcf,['Plots/mean(abs(corr(',GROUP_NAME,'_rcor', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'_grp',...
      num2str(groups(1)),':',num2str(groups(end)),'_calwdw',...
      num2str(calwdw(1)),':',num2str(calwdw(end)),'.jpg'                 ])
  
% More Plotting
for rcorwdw = 1:3
    
subplot(3,2,rcorwdw*2-1)
GROUP_NAME = 'glb_ts';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rgbmap = hsv(7);
for i=1:7
    plot(squeeze(all_expvar_EPC(rcorwdw,i*10,calwdw,:)),'Color',rgbmap(i,:)); hold on;
end
hold off; legend('group=10','group=20','group=30','group=40','group=50','group=60','group=70','location','northeast');
xlabel('EOF'); ylabel('Percentage Variance Explained');
title(['Var Explained ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs)'   ])
  
subplot(3,2,rcorwdw*2)
GROUP_NAME = 'glb_ts_nstat_sigpcd';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rgbmap = hsv(7);
for i=1:7
    plot(squeeze(all_expvar_EPC(rcorwdw,i*10,calwdw,:)),'Color',rgbmap(i,:)); hold on;
end
hold off; legend('group=10','group=20','group=30','group=40','group=50','group=60','group=70','location','northeast');
xlabel('EOF'); ylabel('Percentage Variance Explained');
title(['Var Explained ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs),'  ])
end

set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]);
% saveas(gcf,['Plots/varexp_EOF_',GROUP_NAME,'_rcor', ...
%       num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),...
% '.jpg'                 ])

% More Plotting EOF vs n34
for rcorwdw = 1:3
    
subplot(3,2,rcorwdw*2-1)
GROUP_NAME = 'glb_ts';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rgbmap = hsv(7);
for i=1:7
    plot(squeeze(all_PC_EPC(rcorwdw,i*10,calwdw,:)),'Color',rgbmap(i,:)); hold on;
end
hold off; legend('group=10','group=20','group=30','group=40','group=50','group=60','group=70','location','northeast');
xlabel('EOF'); ylabel('Percentage Variance Explained');
title(['Var Explained ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs)'   ])
  
subplot(3,2,rcorwdw*2)
GROUP_NAME = 'glb_ts_nstat_sigpcd';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rgbmap = hsv(7);
for i=1:7
    plot(squeeze(all_expvar_EPC(rcorwdw,i*10,calwdw,:)),'Color',rgbmap(i,:)); hold on;
end
hold off; legend('group=10','group=20','group=30','group=40','group=50','group=60','group=70','location','northeast');
xlabel('EOF'); ylabel('Percentage Variance Explained');
title(['Var Explained ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs),'  ])
end

set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]);

%% Plotting PC time series
calwdw=5;
for rcorwdw = 1:3
    
subplot(3,1,rcorwdw)
GROUP_NAME = 'glb_ts_nstat_sigpcd';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rgbmap = hsv(7); corr_PC_n34 = nan(1,7);
for i=1:7
    plot(squeeze(all_PC_CPS(rcorwdw,i*10,calwdw,1,:)),'Color',rgbmap(i,:)); hold on;
    corr_PC_n34(i) = corr(squeeze(all_PC_EPC(rcorwdw,i*10,calwdw,1,:)),n34_ind);
end
hold off; legend(['g=10, cor=',num2str(corr_PC_n34(1),'%.2f')], ...
                 ['g=20, cor=',num2str(corr_PC_n34(2),'%.2f')], ...
                 ['g=30, cor=',num2str(corr_PC_n34(3),'%.2f')], ...
                 ['g=40, cor=',num2str(corr_PC_n34(4),'%.2f')], ...
                 ['g=50, cor=',num2str(corr_PC_n34(5),'%.2f')], ...
                 ['g=60, cor=',num2str(corr_PC_n34(6),'%.2f')], ...
                 ['g=70, cor=',num2str(corr_PC_n34(7),'%.2f')], ...
                 'location','northeast');
xlabel('Year');
title(['PC1 Time series ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw(1)),':',num2str(rcorwdw(end)),'(',num2str(rcorwdw*30+1),'yrs), calwdw=', num2str(calwdw) ...
      ])
end
set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]);

%% Regression of EOF1 to ats
GROUP_NAME = 'glb_ts';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
rcorwdw=3; groups=68; calwdw=9; 
ats_reg_eof1 = nan(2,size(ats,2),size(ats,3));
if corr(squeeze(all_PC_EPC(rcorwdw,groups,calwdw,1,:)),n34_ind) < 0
    PC_temp = -squeeze(all_PC_EPC(rcorwdw,groups,calwdw,1,:));
else
    PC_temp = squeeze(all_PC_EPC(rcorwdw,groups,calwdw,1,:));
end
    
for i=1:length(lat)
    for j=1:length(lon)
        ats_reg_eof1(:,i,j) = regress(squeeze(ats(:,i,j)), [PC_temp,ones(499,1)]);
    end
end
pcolor(lon,lat,squeeze(ats_reg_eof1(1,:,:))); plotworld;
colormap(b2r(-0.05,0.05)); colorbar;
title(['Reg coefs - reconstruction PC1 and model temp ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw),'(',num2str(rcorwdw*30+1),'yrs), group=',...
      num2str(groups),', calwdw=',...
      num2str(calwdw)                                ])
  
%% Plotting Iterations vs Groups for EOF1 values

calwdw=10;
GROUP_NAME = 'glb_ts_nstat_sigpcd';
load(['DataFiles/reconEOFs_',GROUP_NAME,'.mat']);
for i=1:3
    subplot(3,1,i)
pcolor(double(abs(squeeze(all_eof_EPC(i,:,calwdw,1,:))))); shading flat;
colorbar; colormap(flipud(hot)); caxis([0 0.03]);
title(['EOF1 Values ',strrep(GROUP_NAME,'_','\_'),': rcorwdw=', ...
      num2str(rcorwdw),'(',num2str(rcorwdw*30+1),'yrs), calwdw=',...
      num2str(calwdw)                                ])
xlabel('Reconstruction Iteration'); ylabel('Group no. (no. of proxies in reconstruction)');
end

%% Finding bad iterations

all_eof1_EPC = squeeze(all_eof_EPC(:,:,:,1,:));
[~,lowest_EOF] = min(abs(all_eof1_EPC(:)));
[a,b,c,d] = ind2sub([3,70,10,1000],lowest_EOF)
corr(all_stn_EPC(a,b,c,1,d),n34_ind)

%% 