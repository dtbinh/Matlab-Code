% This script plots a quantile quantile plot of the model correlations over
% 500 years, vs the running correlation for the same period.

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
corr_pr = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_pr(i,j) = corr(n34_ind,apr(:,i,j));
    end
end

clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

%% Loading Stuff
window = 31; NUM_YRS = 499;
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

%% Line Quantile Plots (newer version)

% Temp
bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_ts(:));
sorted_ts_runcorr = ts_runcorr(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; ts_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    ts_runcorr_quan(:,m) = quantile(reshape(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% The plotting part

clf; axes; hold on; Hnd = nan(1,7); plot([0,-1;0,1],[-1,0;1,0],'k')
for n=1:7
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(text(0,0.30,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.07,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',30,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['Running Correlations (',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled temperature, rcor=',num2str(window),'yrs'])
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/quant_plots_ts_rcor',num2str(window),'yr.jpg'])

% Precip

bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_pr(:));
sorted_pr_runcorr = pr_runcorr(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; pr_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    pr_runcorr_quan(:,m) = quantile(reshape(sorted_pr_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% The plotting part

clf; axes; hold on; Hnd = nan(1,7); plot([0,-1;0,1],[-1,0;1,0],'k')
for n=1:7
    HA(n) = plot(bin,pr_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(HA([3,5]),'Visible','off')
set(text(0,0.31,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0.0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.06,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',30,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.35,-0.40,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['Running Correlations (',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled precipitation, rcor=',num2str(window),'yrs'])
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/quant_plots_pr_rcor',num2str(window),'yr.jpg'])
%% Plot of Variability in Correlation
% Rainfall
clf;
hold on;
for i=1:size(pr_quan,1)
    scatter(corr_pr(:),pr_quan(i,:),1,nonstat_prmap(:));
    colormap(flipud(hot))
    caxis([0 200])
end
hold off;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
axis equal
axis([-1,1,-1,1]);
grid on
xlabel('Mean correlation for entire dataset');
ylabel('30 year running correlation coefficients');
title('Spread of correlation quantiles for Precipitation');
h=colorbar;
title(h,'No. of Nonstationary Years')

% Temperature
figure;
clf;
hold on;
for i=1:size(ts_quan,1)
    scatter(corr_ts(:),ts_quan(i,:),1,nonstat_tsmap(:));
    colormap(flipud(hot))
    caxis([0 200])
end
hold off;
% h=legend('0','0.05','0.25','0.5','0.75','0.95','1','Location','SouthEast');
% title(h,'Percentiles');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
axis equal
axis([-1,1,-1,1]);
grid on
xlabel('Mean correlation for entire dataset');
ylabel('30 year running correlation coefficients');
title('Spread of correlation quantiles for Temperature');
h=colorbar;
title(h,'No. of Nonstationary Years')

%% Plot of nonstat dots

clf; axes; hold on;
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
for l=(length(ts_runcorr(:,1,1))+1-30-31):-1:1 % Running nonstat window and runcorr window
    sig_stns = find(running_nonstat_tsmap(l,:)>=10);
    scatter(corr_ts(sig_stns),squeeze(ts_runcorr(31+l+15,sig_stns)),1,squeeze(running_nonstat_tsmap(l,sig_stns)))
    colormap(flipud(gray));
    ylabel('Running Correlations (30 yr windows)'); xlabel('Correlations over 499 yr period');
    title('Correlation percentiles for modeled temperature')
    grid on; axis equal; xlim([-1 1]); ylim([-1 1]);
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
colorbar; caxis([0 30]);
saveas(gcf,'Plots/scatter(corr_ts,ts_runcorr,run_nstat_map).jpg')

figure; clf; axes; hold on;
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
for l=(length(pr_runcorr(:,1,1))+1-30-31):-1:1 % Running nonstat window and runcorr window
    sig_stns = find(running_nonstat_prmap(l,:)>=10);
    scatter(corr_pr(sig_stns),squeeze(pr_runcorr(31+l+15,sig_stns)),1,squeeze(running_nonstat_prmap(l,sig_stns)))
    colormap(flipud(gray));
    ylabel('Running Correlations (30 yr windows)'); xlabel('Correlations over 499 yr period');
    title('Correlation percentiles for modeled precipitation')
    grid on; axis equal; xlim([-1 1]); ylim([-1 1]);
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
colorbar; caxis([0 30]);
saveas(gcf,'Plots/scatter(corr_pr,pr_runcorr,run_nstat_map).jpg')
%% Plot of non-running nonstat dots,

figure; clf; axes; hold on;
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
sig_stns = find(nonstat_tsmaprecord>0);
corr_3d = permute(repmat(corr_ts,[1 1 499]),[3 1 2]);
scatter(corr_3d(sig_stns),ts_runcorr(sig_stns),1,squeeze(nonstat_tsmaprecord(sig_stns)),'k')
ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled temperature - rcor:',num2str(window),'yr'])
grid on; axis equal; xlim([-1 1]); ylim([-1 1]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/scatter(corr_ts,ts_runcorr,nstat_map)_rcor',num2str(window),'.jpg'])
close;

figure; clf; axes; hold on;
plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
sig_stns = find(nonstat_prmaprecord>0);
corr_pr_3d = permute(repmat(corr_pr,[1 1 499]),[3 1 2]);
scatter(corr_pr_3d(sig_stns),pr_runcorr(sig_stns),1,squeeze(nonstat_prmaprecord(sig_stns)),'k')
ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled precipitation - rcor:',num2str(window),'yr'])
grid on; axis equal; xlim([-1 1]); ylim([-1 1]);
set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/scatter(corr_pr,pr_runcorr,nstat_map)_rcor',num2str(window),'.jpg'])
close;

%% Plotting Density plot 

for window=[31,61,91]
    load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

    clf; axes; hold on;
    plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
    plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
    sig_stns = find(nonstat_tsmaprecord>0);
    corr_3d = permute(repmat(corr_ts,[1 1 499]),[3 1 2]);
    data1 = corr_3d(sig_stns);
    data2 = ts_runcorr(sig_stns);
    values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1}); % dont know what the 51s are
    colormap(flipud(gray))
    caxis([0,200])
    % imagesc(values)
    pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
    colorbar
    axis equal
    xlim([-1,1]); ylim([-1, 1]);
    ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
    title(['Correlation percentiles for modeled temperature - rcor:',num2str(window),'yr'])
    set(gcf, 'PaperUnits', 'centimeters'); % May already be default
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    saveas(gcf,['Plots/hist3(corr_ts_3d,ts_runcorr)_rcor',num2str(window),'.jpg'])
    xlim([-1,1]); ylim([-1, 1]); hold off;

    clf; axes; hold on;
    plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
    plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
    sig_stns = find(nonstat_prmaprecord>0);
    corr_pr_3d = permute(repmat(corr_pr,[1 1 499]),[3 1 2]);
    data1 = corr_pr_3d(sig_stns);
    data2 = pr_runcorr(sig_stns); 
    values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1}); % dont know what the 51s are
    colormap(flipud(gray))
    caxis([0,200])
    % imagesc(values)
    pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
    colorbar
    axis equal
    xlim([-1,1]); ylim([-1, 1]);
    ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
    title(['Correlation percentiles for modeled precipitation - rcor:',num2str(window),'yr'])
    set(gcf, 'PaperUnits', 'centimeters'); % May already be default
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    saveas(gcf,['Plots/hist3(corr_pr_3d,pr_runcorr)_rcor',num2str(window),'.jpg'])
    xlim([-1,1]); ylim([-1, 1]); hold off;
    close;
    
end

%% Plotting Density plot with other conditions

for window=[31,61,91]
    load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])

    clf; axes; hold on;
    plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
    plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
    
    nonstat_map_ind = find(nonstat_tsmap > ceil(0.1*(NUM_YRS-window)) & ...
        1 );
    % squeeze(abs(mean(ts_pc(1,(window+2):end,:,:)-ts_pc(2,(window+2):end,:,:),2))) > 0.3
    [a,b] = ind2sub(size(nonstat_tsmap),nonstat_map_ind);
    corr_3d = permute(repmat(corr_ts,[1 1 499]),[3,1,2]);
    runcr = nan(size(ts_runcorr)); corr_3d_fmt = nan(size(corr_3d));
    for i=1:length(a)
        indices = find(nonstat_tsmaprecord(:,a(i),b(i)));
        corr_3d_fmt(indices,a(i),b(i)) = corr_3d(indices,a(i),b(i));
        runcr(indices,a(i),b(i)) = ts_runcorr(indices,a(i),b(i));
    end
    data1 = corr_3d_fmt;
    data2 = runcr;
    values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1});
    colormap(flipud(gray))
    caxis([0,100])
    pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
    colorbar
    axis equal
    xlim([-1,1]); ylim([-1, 1]);
    ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
    title(['Correlation percentiles for modeled temperature - rcor:',num2str(window),'yr'])
    set(gcf, 'PaperUnits', 'centimeters'); % May already be default
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    saveas(gcf,['Plots/hist3(corr_ts_3d,ts_runcorr)_rcor_nstat',num2str(window),'.jpg'])
    xlim([-1,1]); ylim([-1, 1]); hold off;

    clf; axes; hold on;
    plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
    plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
    nonstat_map_ind = find(nonstat_prmap > ceil(0.1*(NUM_YRS-window)) & ...
        1 );
    % squeeze(abs(mean(pr_pc(1,(window+2):end,:,:)-pr_pc(2,(window+2):end,:,:),2))) > 0.3
    [a,b] = ind2sub(size(nonstat_prmap),nonstat_map_ind);
    corr_3d = permute(repmat(corr_pr,[1 1 499]),[3,1,2]);
    runcr = nan(size(pr_runcorr)); corr_3d_fmt = nan(size(corr_3d));
    for i=1:length(a)
        indices = find(nonstat_prmaprecord(:,a(i),b(i)));
        corr_3d_fmt(indices,a(i),b(i)) = corr_3d(indices,a(i),b(i));
        runcr(indices,a(i),b(i)) = pr_runcorr(indices,a(i),b(i));
    end
    data1 = corr_3d_fmt;
    data2 = runcr;
    values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1});
    colormap(flipud(gray))
    caxis([0,100])
    pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
    colorbar
    axis equal
    xlim([-1,1]); ylim([-1, 1]);
    ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
    title(['Correlation percentiles for modeled precipitation - rcor:',num2str(window),'yr'])
    set(gcf, 'PaperUnits', 'centimeters'); % May already be default
    set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
    saveas(gcf,['Plots/hist3(corr_pr_3d,pr_runcorr)_rcor_nstat',num2str(window),'.jpg'])
    xlim([-1,1]); ylim([-1, 1]); hold off;
    close;
    
end

%% Density plots overlayed with Line correlation percentile plots
for window=[31,61,91]
    load(['DataFiles/runcorr',num2str(window),'yrwdw.mat'])
    load(['DataFiles/nonstat_map',num2str(window),'yrwdw.mat'])
% Temp
bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_ts(:));
sorted_ts_runcorr = ts_runcorr(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; ts_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    ts_runcorr_quan(:,m) = quantile(reshape(sorted_ts_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

clf; axes; hold on; Hnd = nan(1,7); plot([0,-1;0,1],[-1,0;1,0],'k')
for n=1:7
    HA(n) = plot(bin,ts_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(HA(4),'LineWidth',3); 

plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
nonstat_map_ind = find(nonstat_tsmap > ceil(0.1*(NUM_YRS-window)) & ...
    1 );
% squeeze(abs(mean(ts_pc(1,(window+2):end,:,:)-ts_pc(2,(window+2):end,:,:),2))) > 0.3
[a,b] = ind2sub(size(nonstat_tsmap),nonstat_map_ind);
corr_3d = permute(repmat(corr_ts,[1 1 499]),[3,1,2]);
runcr = nan(size(ts_runcorr)); corr_3d_fmt = nan(size(corr_3d));
for i=1:length(a)
    indices = find(nonstat_tsmaprecord(:,a(i),b(i)));
    corr_3d_fmt(indices,a(i),b(i)) = corr_3d(indices,a(i),b(i));
    runcr(indices,a(i),b(i)) = ts_runcorr(indices,a(i),b(i));
end
data1 = corr_3d_fmt;
data2 = runcr;
values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1});
colormap(flipud(gray))
pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
colorbar
caxis([0,100])
axis equal
xlim([-1,1]); ylim([-1, 1]);
ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled temperature - rcor:',num2str(window),'yr'])
set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/hist3(corr_ts_3d,ts_runcorr)_rcor_nstat+lines',num2str(window),'.jpg'])
xlim([-1,1]); ylim([-1, 1]); hold off;




% Precip

bin = -1.0:0.02:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_pr(:));
sorted_pr_runcorr = pr_runcorr(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; pr_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    pr_runcorr_quan(:,m) = quantile(reshape(sorted_pr_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% The plotting part

clf; axes; hold on; Hnd = nan(1,7); plot([0,-1;0,1],[-1,0;1,0],'k')
for n=1:7
    HA(n) = plot(bin,pr_runcorr_quan(n,:));
    set(HA(n),'Color','k','LineWidth',2);
end
grid on; axis equal; axis([-1 1 -1 1]);
set(HA([1,7]),'LineStyle','-','LineWidth',1);
set(HA([3,5]),'Visible','off')
set(HA(4),'LineWidth',3);

plot([-0.3,-0.3],[-1 1],'k'); plot([0.3,0.3],[-1 1],'k')
plot([-1 1],[-0.3,-0.3],'k'); plot([-1 1],[0.3,0.3],'k')
nonstat_map_ind = find(nonstat_prmap > ceil(0.1*(NUM_YRS-window)) & ...
    1 );
% squeeze(abs(mean(pr_pc(1,(window+2):end,:,:)-pr_pc(2,(window+2):end,:,:),2))) > 0.3
[a,b] = ind2sub(size(nonstat_prmap),nonstat_map_ind);
corr_3d = permute(repmat(corr_pr,[1 1 499]),[3,1,2]);
runcr = nan(size(pr_runcorr)); corr_3d_fmt = nan(size(corr_3d));
for i=1:length(a)
    indices = find(nonstat_prmaprecord(:,a(i),b(i)));
    corr_3d_fmt(indices,a(i),b(i)) = corr_3d(indices,a(i),b(i));
    runcr(indices,a(i),b(i)) = pr_runcorr(indices,a(i),b(i));
end
data1 = corr_3d_fmt;
data2 = runcr;
values = hist3([data1(:) data2(:)],{-1:0.01:1, -1:0.01:1});
colormap(flipud(gray))
caxis([0,100])
pcolor(-1:0.01:1,-1:0.01:1,values'); shading flat
colorbar
axis equal
xlim([-1,1]); ylim([-1, 1]);
ylabel(['Running Correlations (using ',num2str(window),' yr windows)']); xlabel('Correlations over 499 yr period');
title(['Correlation percentiles for modeled precipitation - rcor:',num2str(window),'yr'])
set(gcf, 'PaperUnits', 'centimeters'); % May already be default
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/hist3(corr_pr_3d,pr_runcorr)_rcor_nstat+lines',num2str(window),'.jpg'])
xlim([-1,1]); ylim([-1, 1]); hold off;
close;
end