%%                      INFORMATION

%  Author: Laur Ferris (lnferris@alum.mit.edu)

% Dependencies: readcacheslocum.m, readslocumbd.m (Robert Todd)

% Engineering Outputs: dbdtime, dbdlon, dbdlat 
% Science Outputs: tt, llon, llat, pprs, ttem, ccon, ddep, SP, SA, CT, rho, sigma0

% The dbd paths can be DBD or SBD. The ebd paths can be EBD or TBD.

%%                      USER INPUTS

%dbd_path = '/Users/gong/Box/glider/gliderData/amelia-20200825-hurricane-realtime/amelia-from-glider-20200911T035349/*.sbd';
%ebd_path = '/Users/gong/Box/glider/gliderData/amelia-20200825-hurricane-realtime/amelia-from-glider-20200911T035349/*.tbd';

dbd_path = '/Users/gong/Box/glider/gliderData/sylvia-20200903-hurricane-realtime/SBD/*.sbd';
ebd_path = '/Users/gong/Box/glider/gliderData/sylvia-20200903-hurricane-realtime/TBD/*.tbd';

% Specify cache file listed in first DBD file (e.g. 'sensor_list_crc: 574a5509').
%cache_path_dbd = '/Users/gong/Box/glider/gliderData/amelia-20200825-hurricane-realtime/cache/be93efad.cac';
cache_path_dbd = '/Users/gong/GitHub/data_processing/cache/be93efad.cac';

% Specify cache file listed in first EBD file (e.g. 'sensor_list_crc: 357d1ddc').
%cache_path_ebd =  '/Users/gong/Box/glider/gliderData/amelia-20200825-hurricane-realtime/cache/488fee97.cac';
cache_path_ebd =  '/Users/gong/GitHub/data_processing/cache/2804f449.cac';

%%                      SCRIPT

% Load raw DBD data.
[dbdtime,gpslat,gpslon] = loadDBD(dbd_path,cache_path_dbd);

% Load raw EBD data.
[ebdtime,con,prs,tem] = loadEBD(ebd_path,cache_path_ebd);

% Process raw GPS.
[gpslat,gpslon] = processGPS(gpslat,gpslon);

% Sort raw data by time.
[dbdtime,ebdtime,gpslat,gpslon,con,prs,tem] = timeSort(dbdtime,ebdtime,gpslat,gpslon,con,prs,tem);

% Spatially interpolate GPS data.
[dbdlat,dbdlon,ebdlat,ebdlon] = interpGPS(dbdtime,ebdtime,gpslat,gpslon);

% Process raw CTD data.
[con,prs,tem] = processCTD(con,prs,tem);

% Practical salinity from conductivity.
sal = gsw_SP_from_C(con,tem,prs);

%% Plot.

figure
subplot(1,2,1)
scatter(ebdtime,-prs,[],tem,'filled')
colorbar
datetick
cmocean('thermal')

subplot(1,2,2)
scatter(ebdtime,-prs,[],sal,'filled')
colorbar
datetick
cmocean('haline')

figure
scatter(ebdlon,ebdlat,'filled')
grid on


%%                      FUNCTIONS

function [dbdtime,gpslat,gpslon] = loadDBD(dbd_path,cache_path_dbd)
% Initialize variables to hold all DBD data in the path.
dbdtime = [];
%z = []; 
%leak = []; 
%vac = []; 
%pitch = []; 
%lon = []; 
%lat = []; 
%roll = []; 
%bat = [];
gpslat = []; 
gpslon = []; 
%yrr = []; 
%mon = []; 
%day = []; 
%hrr = []; 
%mnn = []; 
%sec = [];

% For each DBD file...
data_glob = dir(dbd_path);
cache_glob = dir(cache_path_dbd);
sensor_label = cache_glob(1).name(1:end-4);
sensor_list = readcacheslocum([cache_glob(1).folder '/' cache_glob(1).name]);
cache_directory = [cache_glob(1).folder '/' ];

for i = 1:length(data_glob)
    
    % Build a data structure from the file.
    filename = [data_glob(i).folder '/' data_glob(i).name];  
    data_struct = readslocumbd(filename,sensor_label,sensor_list,cache_directory);
      
    % Write the data structure into previously-initilaized variables.
    disp(['Processing DBD file ', num2str(i), ' of ', num2str(length(data_glob))]);

    dbdtime = [dbdtime; datenum(0,0,0,0,0,data_struct.m_present_time) + datenum(1970,1,1)];
    %z = [z; data_struct.m_depth];
    %leak = [leak; data_struct.m_leakdetect_voltage];
    %vac = [vac; data_struct.m_vacuum];
    %lon = [lon; data_struct.m_lon];
    %lat = [lat; data_struct.m_lat];
    gpslat = [gpslat; data_struct.m_gps_lat];
    gpslon = [gpslon; data_struct.m_gps_lon];
    %pitch = [pitch; data_struct.m_pitch];
    %roll = [roll; data_struct.m_roll];
    %bat = [bat; data_struct.m_battpos];
    %yrr = [yrr; data_struct.m_gps_utc_year];
    %mon = [mon; data_struct.m_gps_utc_month];
    %day = [day; data_struct.m_gps_utc_day];
    %hrr = [hrr; data_struct.m_gps_utc_hour];
    %mnn = [mnn; data_struct.m_gps_utc_minute];
    %sec = [sec; data_struct.m_gps_utc_second];
        
end
disp('Finished processing DBD files.')

end

function [ebdtime,con,prs,tem] = loadEBD(ebd_path,cache_path_ebd)
                     
% Initialize variables to hold all EBD data in the path.
ebdtime = [];
prs = []; 
tem = []; 
con = [];

% For each EBD file...
data_glob = dir(ebd_path);
cache_glob = dir(cache_path_ebd);
sensor_label = cache_glob(1).name(1:end-4);
sensor_list = readcacheslocum([cache_glob(1).folder '/' cache_glob(1).name]);
cache_directory = [cache_glob(1).folder '/' ];

for i = 1:length(data_glob)
    try
    
    % Build a data structure from the file.
    filename = [data_glob(i).folder '/' data_glob(i).name];  
    data_struct = readslocumbd(filename,sensor_label,sensor_list,cache_directory);
      
    % Write the data structure into previously-initilaized variables.
    disp(['Processing EBD file ', num2str(i), ' of ', num2str(length(data_glob))]);
    
    % Catch and ignore files with time-CTD length mismatches.
    if length(data_struct.sci_m_present_time) == length(data_struct.sci_water_pressure)  
        ebdtime = [ebdtime; datenum(0,0,0,0,0,data_struct.sci_m_present_time) + datenum(1970,1,1)];
        prs=[prs;data_struct.sci_water_pressure];
        tem=[tem;data_struct.sci_water_temp];
        con=[con;data_struct.sci_water_cond];
    else 
        disp([i,length(data_struct.sci_m_present_time),length(data_struct.sci_water_pressure)] ) 
        padding = NaN(length(data_struct.sci_m_present_time)-length(data_struct.sci_water_pressure),1); 
        ebdtime = [ebdtime; datenum(0,0,0,0,0,data_struct.sci_m_present_time) + datenum(1970,1,1)];
        prs=[prs;data_struct.sci_water_pressure;padding];
        tem=[tem;data_struct.sci_water_temp;padding];
        con=[con;data_struct.sci_water_cond;padding];     
    end
    end
end
disp('Finished processing EBD files.')
end        

function [gpslat,gpslon] = processGPS(gpslat,gpslon)
% Convert to decimal degrees and remove impossible values.
lat_sign = sign(gpslat);
lon_sign = sign(gpslon);
gpslat = abs(gpslat);
gpslon = abs(gpslon);
gpslat = (floor(gpslat/100)+100*((gpslat/100)-floor(gpslat/100))/60).*lat_sign;
gpslon = (floor(gpslon/100)+100*((gpslon/100)-floor(gpslon/100))/60).*lon_sign;
gpslat(abs(gpslat)>90) = NaN;
gpslon(abs(gpslon)>360) = NaN;
disp('GPS processing applied.')
end

function [dbdtime,ebdtime,gpslat,gpslon,con,prs,tem] = timeSort(dbdtime,ebdtime,gpslat,gpslon,con,prs,tem)
% Sort relevant DBD data by time.
[dbdtime,dbdtime_inds] = sort(dbdtime);
gpslat = gpslat(dbdtime_inds);
gpslon = gpslon(dbdtime_inds);

% Sort relevant EBD data by time.
[ebdtime,ebdtime_inds] = sort(ebdtime);
con = con(ebdtime_inds);
prs = prs(ebdtime_inds);
tem = tem(ebdtime_inds);
end

function [dbdlat,dbdlon,ebdlat,ebdlon] = interpGPS(dbdtime,ebdtime,gpslat,gpslon)

% Find non-NaN GPS indicies for each unique timestep.
[~,unique_dbdtime_inds] = unique(dbdtime);
good_inds = intersect(find(~isnan(gpslon)),unique_dbdtime_inds);

% Remap GPS fixes into complex plane.
gps_zll = gpslon(good_inds) + gpslat(good_inds)*1j; % z = x+iy

% Interpolate good lon/lat data to full data timeseries.
%   vq = interp1(x,v,xq)
%   x is sample points, given by good GPS fixes.
%   v are corresponding values, given by complex lon/lat.
%   xq are query points, given by the DBD or EBD time series.
dbdlonlat = interp1(dbdtime(good_inds),gps_zll,dbdtime); 
ebdlonlat = interp1(dbdtime(good_inds),gps_zll,ebdtime); 

% Bring GPS fixes back into cartesian plane.
dbdlon = real(dbdlonlat); 
dbdlat = imag(dbdlonlat);
ebdlon = real(ebdlonlat); 
ebdlat = imag(ebdlonlat); 
clear gpslat gpslon

disp('GPS interpolation applied.')
end

function [con,prs,tem] = processCTD(con,prs,tem)
% Basic processing.
prs = prs*10; % bars -> decibars
con = con*10; % S/m -> mS/cm.
badind = find(con < 0.1 | tem == 0.0000000);
con(badind) = NaN;
tem(badind) = NaN;
prs(badind) = NaN;
disp('CTD corrections applied.')
end

function [tt,llat,llon,ccon,pprs,ttem] = timeInterp(ebdtime,ebdlat,ebdlon,con,tem,prs)

% Create a fixed 0.5-sec time grid.
interval = datenum(0,0,0,0,0,0.5); 
tt = ceil(ebdtime(find(~isnan(ebdtime),1,'first'))):interval:floor(ebdtime(find(~isnan(ebdtime),1,'last')));

% Interpolate the data to the fixed time grid, ommitting > 5-sec gaps for CTD data.
t_gap = datenum(0,0,0,0,0,5);
llon = interp1(ebdtime(~isnan(ebdlon)),ebdlon(~isnan(ebdlon)),tt,'linear');
llat = interp1(ebdtime(~isnan(ebdlat)),ebdlat(~isnan(ebdlat)),tt,'linear');
pprs = interp1gap(ebdtime(~isnan(prs)),prs(~isnan(prs)),tt,t_gap,'linear');
ttem = interp1gap(ebdtime(~isnan(tem)),tem(~isnan(tem)),tt,t_gap,'linear');
ccon = interp1gap(ebdtime(~isnan(con)),con(~isnan(con)),tt,t_gap,'linear');

disp('Time interpolation applied.')
end


 function [ddep,SP,SA,CT,rho,sigma0] = derivedParams(llat,llon,ccon,pprs,ttem)
% Depth [m] below surface.
ddep = gsw_z_from_p(pprs,llat);

% Practical salinity from conductivity.
SP = gsw_SP_from_C(ccon,ttem,pprs);

% Absolute salinity from practical salinity.
[SA, ~] = gsw_SA_from_SP(SP,pprs,llon,llat);

% Conservative temperature.
CT = gsw_CT_from_t(SA,ttem,pprs);

% In situ density.
rho = gsw_rho(SA,CT,pprs);

% Potential density anomaly.
sigma0 = gsw_sigma0(SA,CT);

disp('Derived parameters calculated.')
end

