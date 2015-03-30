%% This script calculates the Nino3.4 teleconnection strengths on GISTEMP1200_ERSST instrumental data

%% Date Conversions
START_YEAR = 1800;
START_MONTH = 1;
START_DAY = 1;
FILE_NAME = 'DataFiles/gistemp1200_ERSST.nc';
DATE_OFFSET = datenum(START_YEAR,START_MONTH,START_DAY);
sst_datenum = nc_varget(FILE_NAME,'time')+DATE_OFFSET;

%% Nino 3.4 SSTs
ssta = nc_varget(FILE_NAME,'tempanomaly');
lat = nc_varget(FILE_NAME,'lat');
lon = nc_varget(FILE_NAME,'lon');
ssta_datenum = nc_varget(FILE_NAME,'time')+DATE_OFFSET;

% Changing SSTAs to 0-360E
a=ssta(:,:,1:90);
ssta(:,:,1:90)=ssta(:,:,91:end);
ssta(:,:,91:end)=a;
lon=lon+180;

% Nino 3.4 Indexes
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nW] = min(abs(lon-190));
[~,nE] = min(abs(lon-240));

% Annual (Jul-Jun) Means and Anomalies
% jul_jun_fmt = 7:1614;
% ssta=ssta(jul_jun_fmt,:,:);
% time=ssta_datenum(jul_jun_fmt);

windowSize = 12;
F = ones(1,windowSize)/windowSize;
trend = zeros(size(ssta,1),size(ssta,2),size(ssta,3));
for i=1:size(ssta,2)
    for j=1:size(ssta,3)
        trend(:,i,j) = filter(F,1,squeeze(ssta(:,i,j)));
        ssta(7:end,i,j) = ssta(7:end,i,j) - trend(1:end-6,i,j);
    end
end

% Annual (Jul-Jun) Means and Anomalies
jul_jun_fmt = 7:1614;
ssta=ssta(jul_jun_fmt,:,:);
time=ssta_datenum(jul_jun_fmt);


ssta_an=zeros(size(ssta,1)/12,size(ssta,2),size(ssta,3));
for y=1:length(time)/12
    ssta_an(y,:,:)=mean(ssta((12*(y-1)+1):(y*12),:,:),1);
end

n34_ind = mean(mean(ssta_an(:,nS:nN,nW:nE),3),2);


%% Calculating Correlations 

corr_ssta_an = zeros(size(ssta,2),size(ssta,3));
for i=1:size(ssta,2)
    for j=1:size(ssta,3)
        corr_ssta_an(i,j) = corr(n34_ind(end-50:end),ssta_an(end-50:end,i,j));
    end
end


contourf(lon,lat,corr_ssta_an); axis equal;
shading flat

title('r(sta,Nino3.4)','FontSize',14,'FontWeight','bold');
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
% hold on; pcolor(lon,lat,double(isnan(corr_ssta_an))); hold off
% print -painters -dpdf -r600 test.pdf
set(gcf, 'PaperSize', [35 13]);
set(gca, 'color', [0 0 0]);
