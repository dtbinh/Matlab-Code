% This script appends the matlab files that contain the synthetic running
% correlations for backup purposes or other reasons. 

%% Setup

ts_file = 'DataFiles/ts_A1.nc'; pr_file = 'DataFiles/pr_A1.nc';

lat = nc_varget(ts_file,'lat');
lon = nc_varget(ts_file,'lon');

YEARS = 499;
LAT_LEN = 90;
LON_LEN = 144;


%% Combining the existing files

for n=1:1000

    temp_prsynruncorr = nan(YEARS, LAT_LEN, LON_LEN);
    temp_tssynruncorr = nan(YEARS, LAT_LEN, LON_LEN);
    
    load(['Synth_corr_31yrwin/run',num2str(n),'syncorr.mat']);
    temp_prsynruncorr = pr_synruncorr;
    temp_tssynruncorr = ts_synruncorr;
    
    % Joining at 120E Boundary
    load(['Synth_corr_31yrwin_0-120E/run',num2str(n),'syncorr.mat']);
    S_lat = -60; N_lat = 60; W_lon = 0; E_lon = 120;
    [~,S_bound]= min(abs(lat-S_lat));
    [~,N_bound]= min(abs(lat-N_lat));
    [~,W_bound]= min(abs(lon-W_lon));
    [~,E_bound]= min(abs(lon-E_lon));
    
    if  ~(abs(temp_prsynruncorr(:,S_bound:N_bound,E_bound) - pr_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001 | ...
            abs(temp_tssynruncorr(:,S_bound:N_bound,E_bound) - ts_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001)
        error('Values do not match at 120E Boundary');
    end
        
    temp_prsynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        pr_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    temp_tssynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        ts_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    
    % Joining at 270E Boundary
    load(['Synth_corr_31yrwin_270-360E/run',num2str(n),'syncorr.mat']);
    S_lat = -60; N_lat = 60; W_lon = 270; E_lon = 360;
    [~,S_bound]= min(abs(lat-S_lat));
    [~,N_bound]= min(abs(lat-N_lat));
    [~,W_bound]= min(abs(lon-W_lon));
    [~,E_bound]= min(abs(lon-E_lon));
    
    if  ~(abs(temp_prsynruncorr(:,S_bound:N_bound,W_bound) - pr_synruncorr(:,S_bound:N_bound,W_bound)) <0.00001 | ...
            abs(temp_tssynruncorr(:,S_bound:N_bound,W_bound) - ts_synruncorr(:,S_bound:N_bound,W_bound)) <0.00001)
        error('Values do not match at 270E Boundary');
    end
    
    temp_prsynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        pr_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    temp_tssynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        ts_synruncorr(:,S_bound:N_bound,W_bound:E_bound);    
    
    % Joining at -60N Boundary
    load(['Synth_corr_31yrwin_polarS/run',num2str(n),'syncorr.mat']);
    S_lat = -90; N_lat = -60; W_lon = 0; E_lon = 360;
    [~,S_bound]= min(abs(lat-S_lat));
    [~,N_bound]= min(abs(lat-N_lat));
    [~,W_bound]= min(abs(lon-W_lon));
    [~,E_bound]= min(abs(lon-E_lon));
    
    if  ~(abs(temp_prsynruncorr(:,S_bound:N_bound,E_bound) - pr_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001 | ...
            abs(temp_tssynruncorr(:,S_bound:N_bound,E_bound) - ts_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001)
        error('Values do not match at -60N Boundary');
    end
        
    temp_prsynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        pr_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    temp_tssynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        ts_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    
    % Joining at 60N Boundary
    load(['Synth_corr_31yrwin_polarN/run',num2str(n),'syncorr.mat']);
    S_lat = 60; N_lat = 90; W_lon = 0; E_lon = 360;
    [~,S_bound]= min(abs(lat-S_lat));
    [~,N_bound]= min(abs(lat-N_lat));
    [~,W_bound]= min(abs(lon-W_lon));
    [~,E_bound]= min(abs(lon-E_lon));
    
    if  ~(abs(temp_prsynruncorr(:,S_bound:N_bound,E_bound) - pr_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001 | ...
            abs(temp_tssynruncorr(:,S_bound:N_bound,E_bound) - ts_synruncorr(:,S_bound:N_bound,E_bound)) <0.00001)
        error('Values do not match at 60N Boundary');
    end
        
    temp_prsynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        pr_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    temp_tssynruncorr(:,S_bound:N_bound,W_bound:E_bound) = ...
        ts_synruncorr(:,S_bound:N_bound,W_bound:E_bound);
    
    pr_synruncorr = temp_prsynruncorr;
    ts_synruncorr = temp_tssynruncorr;
    save(['Synth_corr_all/run',num2str(n),'syncorr.mat'],'pr_synruncorr','ts_synruncorr');

end

 