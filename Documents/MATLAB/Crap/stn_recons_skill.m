
% This script will plot the time series and statisitics of the
% reconstruction of n34 from extrapolated 50 year windows of temperature regressions to
% n34, and then calculate the correlation coefficients. This will be done
% over the region 120-280E, -10 to 10 N.

%% Setup
ts_file = 'ts_A1.nc'; pr_file = 'pr_A1.nc';

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

%% Extrapolation of One Point

x_lon = 150; [~,x_ind]= min(abs(lon-x_lon));
y_lat = 0; [~,y_ind]= min(abs(lat-y_lat));
year_wdw = [1 50; 51 100; 101 150; 151 200; 201 250; 251 300; 301 350; ...
               351 400; 401 450; 451 499];
con_n34 = nan(size(year_wdw,1), 499);

for n=1:size(year_wdw,1)
    
    ts_coef = regress(n34_ind(year_wdw(n,1):year_wdw(n,2)), ...
                      [ats(year_wdw(n,1):year_wdw(n,2), y_ind, x_ind), ...
                       ones(1+year_wdw(n,2)-year_wdw(n,1),1)          ]    );
    con_n34(n,:) = ts_coef(1)*ats(:,y_ind,x_ind)';
    
end

%% Extrapolation of an Area

S_lat = -10; N_lat = 10; W_lon = 120; E_lon = 280;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));




%% Plotting series
for n=1:10
    subplot(5,2,n)
    plot(n34_ind,'k');
    hold on;
    plot(con_n34(n,:))
    hold off;
    title(['At ',num2str(y_lat),'{\circ}N, ',num2str(x_lon'),'{\circ}E', ...
           '  Corr: ' , num2str(corr(n34_ind,con_n34(n,:)'),4), ...
           ', Var diff is ', num2str(var(n34_ind)-var(con_n34(n,:)),4) 
%           ', RMSE is ',num2str(sqrt(mean((n34_ind-con_n34(n,:)').^2)),4), ...
           ]);
    ylim([-3.5,3.5])
end

%% Plotting running Variance of Series

n34_mv = movingvar(n34_ind,MV_WDW);
MV_WDW = 20;

for n=1:10

    con_mv = movingvar(con_n34(n,:)',MV_WDW);
    subplot(5,2,n)
    plot(n34_mv,'k');
    hold on;
    plot(con_mv);
    hold off;
    title(['At ',num2str(y_lat),'{\circ}N, ',num2str(x_lon'),'{\circ}E', ...
           '  Corr: ' , num2str(corr(n34_mv(10:end-10),con_mv(10:end-10)),4), ...
           ', Var diff is ', num2str(var(n34_ind)-var(con_n34(n,:)),4), ...
           ', RMSE is ',num2str(sqrt(mean((n34_mv-con_mv).^2)),4) ...
           ]);
    ylim([0,2])
end

% select 50 year window from temp data, and regress against N34

