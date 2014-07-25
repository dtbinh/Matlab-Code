% This script will examine the non-linearity of a gridpoint's response to
% ENSO

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

clear ats_anmn apr_anmn trend ts pr time jul_jun_fmt nN nE nS nW ts_file pr_file i j y

%% Plotting
% ats_n = (ats(:,1,1)./std(ats(:,1,1)));
% apr_n = (apr(:,1,1)./std(apr(:,1,1)));
% n34_ind_n = std((n34_ind./std(n34_ind)));
[sorted,sort_ind] = sort(squeeze(n34_ind));
[~,y]= min(abs(lat--20)); y=45;
[~,x]= min(abs(lon-276)); x=70;
scatter(n34_ind,squeeze(ats(:,y,x)));
% xlim([-3,3]);
% ylim([-2,5]*10^-5);
% xlabel('Nino3.4 Index'); ylabel('ats');
sbd = nan(size(squeeze(ats(1,:,:))));
% find(n34_ind>0)

for i=1:length(lat)
    for j=1:length(lon)
        sbd(i,j) = sum(abs(detrend(squeeze(ats(sort_ind,i,j)))));
%         apr_1(i,j) = squeeze(apr(sort_ind,i,j))
    end
end

pcolor(lon,lat,sbd); plotworld; 
% colormap(b2r(-1*10^-13,1*10^-13));
colormap(flipud(hot))
colorbar;
% caxis([-1e-10,1e-10])

%% Using difference between El Nino slopes and La Nina slopes

[sorted,sort_ind] = sort(squeeze(n34_ind));
plsd = nan(size(squeeze(ats(1,:,:)))); % Phase linear slope difference
nina_lincoef = nan(2,90,144); nino_lincoef = nan(2,90,144);
for i=1:length(lat)
    for j=1:length(lon)
        a_sorted = squeeze(ats(sort_ind,i,j));
        middle = find(sorted>0,1);
        nina_linco = polyfit(sorted(1:middle-1),a_sorted(1:middle-1),1);
        nino_linco = polyfit(sorted(middle:end),a_sorted(middle:end),1);
        nina_lincoef(:,i,j) = nina_linco; nino_lincoef(:,i,j) = nino_linco;
        plsd(i,j) = nino_linco(1)-nina_linco(1);
    end
end

[c,h] = contourf(lon,lat,plsd,10,'b'); plotworld; colorbar; %colormap(redblue(13));
%clabel(c,h,'Color','b')
title('Difference between linear fit slopes of El Nino and La Nina (nino-nina) (temp)');

subplot(2,1,1)
pcolor(lon,lat,nino_lincoef); plotworld;
title('n34>0 regression coefficient')
subplot(2,1,2)
pcolor(lon,lat,nina_lincoef); plotworld;
title('n34<0 regression coefficient')
colormap(b2r(-1.5,1.5))
% colormap(b2r(-5e-5,5e-5))

%% Making ella temperature series

ella_ts = nan(499,90,144);
for i=1:length(lat)
    for j=1:length(lon)
        el_ind = find(n34_ind>=0);
        la_ind = find(n34_ind<0);
        ella_ts(el_ind,i,j) = n34_ind(el_ind)*nino_lincoef(1,i,j)+nino_lincoef(2,i,j);
        ella_ts(la_ind,i,j) = n34_ind(la_ind)*nina_lincoef(1,i,j)+nina_lincoef(2,i,j);
    end
end
        
save('DataFiles/ella_ts.mat', 'ella_ts')    
        
%% Plotting individual points

[~,y]= min(abs(lat-1));
[~,x]= min(abs(lon-153));

a_sorted = squeeze(ats(sort_ind,y,x));
middle = find(sorted>0,1);
scatter(n34_ind,squeeze(ats(:,y,x)),'g'); hold on;
plot(sorted,squeeze(ats(sort_ind,y,x)),'g')
nina_lincoef = polyfit(sorted(1:middle-1),a_sorted(1:middle-1),1);
nino_lincoef = polyfit(sorted(middle:end),a_sorted(middle:end),1);
plot(linspace(-3,0),nina_lincoef(2)+nina_lincoef(1)*linspace(-3,0)); 
plot(linspace(0,3),nino_lincoef(2)+nino_lincoef(1)*linspace(0,3)); hold off;