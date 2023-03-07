function dat = readslocumbd(fname,sensor_list_crc0,sensorlist0,cachedir)

% Read .*bd binary data file from Webb Slocum glider.
% Inputs are:
%   fname: name of file to read (string)
%   sensor_list_crc0: string label for input sensor list
%   sensorlist0: sensor list structure as produced by readcacheslocum
%   cachedir: directory in which cache (.CAC) file resides
% Final three inputs may be left off if file to be read is known to include
% sensor list (i.e., sensor_list_factored == 0).
%
% Output is a structure containing all data from input file with field
% names based on sensor list.
%
% Robert Todd, 14 May 2014

% parameters
numhead = 14;
strheadlist = {'dbd_label','all_sensors','full_filename','filename_extension','mission_name','fileopen_time','sensor_list_crc'};
ncycmax = 1e4;

% open file
fid = fopen(fname,'r');
if fid<0
    fprintf('File NOT Found: %s\n',fname);
    return
end

% read header
for ihead = 1:numhead
    line = fgetl(fid);
    x = textscan(line,'%s %s'); % first is field name, second is value
    fieldname = x{1}{1}(1:end-1);
    
    if any(strcmp(fieldname,strheadlist))
        dat.hdr.(fieldname) = x{2}{1};
    else
        dat.hdr.(fieldname) = str2double(x{2}{1});
    end
end

% sensor list
nsens = dat.hdr.total_num_sensors;
if dat.hdr.sensor_list_factored == 0 % then sensor list present in this file
    dat.sensorlist.transmitted = false(nsens,1);
    dat.sensorlist.sensnum = uint16(zeros(nsens,1)); % will be 0:nsens-1
    dat.sensorlist.indexnum = int16(zeros(nsens,1)); % could be -1 if not transmitted
    dat.sensorlist.nbytes = uint8(zeros(nsens,1)); % 1 = 1 byte integer, 2 = 2 byte integer, 4 = 4 byte float, 8 = 8 byte double
    dat.sensorlist.name = cell(nsens,1);
    dat.sensorlist.unit = cell(nsens,1);
    for isens = 1:nsens
        line = fgetl(fid);
        x = textscan(line,'%*s %s %f %f %f %s %s');
        if strcmp(x{1}{1},'T')
            dat.sensorlist.transmitted(isens) = true;
        end
        dat.sensorlist.sensnum(isens) = uint16(x{2});
        dat.sensorlist.indexnum(isens) = int16(x{3});
        dat.sensorlist.nbytes(isens) = uint8(x{4});
        dat.sensorlist.name{isens} = x{5}{1};
        dat.sensorlist.unit{isens} = x{6}{1};
    end
    
else % sensor list not included in this file
    if exist('sensorlist0','var') && exist('sensor_list_crc0','var') && strcmp(dat.hdr.sensor_list_crc,sensor_list_crc0)
        % then correct sensor list provided as input, use it
        dat.sensorlist = sensorlist0;
    else % need to read .CAC file
        cachename = [cachedir dat.hdr.sensor_list_crc '.CAC'];
        if exist(cachename,'file')
            dat.sensorlist = readcacheslocum(cachename);
        else
            fprintf('Cache file not found for %s\n',fname)
            return
        end
    end
end

% initialize sensor fields
for isens = 1:nsens
    if dat.sensorlist.nbytes(isens) == 8
        dat.(dat.sensorlist.name{isens}) = nan(ncycmax,1,'double');
    else
        dat.(dat.sensorlist.name{isens}) = nan(ncycmax,1,'single');
    end
end

% binary check values
x = fread(fid,1,'*char');
if ~strcmp(x,'s');
    fprintf('Failed binary check value 1: %s\n',fname);
    return
end
x = fread(fid,1,'*char','b');
if ~strcmp(x,'a')
    fprintf('Failed binary check value 2: %s\n',fname);
    return
end
x = fread(fid,1,'*int16','b');
if x ~= 4660
    fprintf('Failed binary check value 3: %s\n',fname);
    return
end
x = fread(fid,1,'*float32','b');
dx = floor(1e3*(x-123.456));
if dx ~= 0
    fprintf('Failed binary check value 4: %s\n',fname);
    return
end
x = fread(fid,1,'*double','b');
if x ~= 123456789.12345
    fprintf('Failed binary check value 5: %s\n',fname);
    return
end

% data cycles
icycle = 0;
lastvalue = cell(nsens,1);
while 1 % keep reading until end of file tag found
    updated = nan(nsens,1);
    x = fread(fid,1,'*char');
    if strcmp(x,'X')  % end of file found
        break
    elseif strcmp(x,'d'); % data cycle tag
        icycle = icycle+1;
        % state bytes
        x = fread(fid,4*dat.hdr.state_bytes_per_cycle,'ubit2','b');
        for isens = 1:nsens
            jsens = dat.sensorlist.indexnum(isens);
            if jsens ~= -1 % then sensor transmitted
                updated(isens) = x(jsens+1); % 0: not updated, 1: updated with same value, 2: updated with new value
            end
        end
        
        % sensor data
        for isens = 1:nsens
            if updated(isens) == 2 % then expect data
                switch dat.sensorlist.nbytes(isens)
                    case 1
                        x = fread(fid,1,'*int8','b');
                    case 2
                        x = fread(fid,1,'*int16','b');
                    case 4
                        x = fread(fid,1,'*single','b');
                    case 8
                        x = fread(fid,1,'*double','b');
                end
                dat.(dat.sensorlist.name{isens})(icycle) = x;
                lastvalue{isens} = x;
                
            elseif updated(isens) == 1 % assign previous value
                dat.(dat.sensorlist.name{isens})(icycle) = lastvalue{isens};
            end
        end
    else
        fprintf('Problem reading binary data: %s\n',fname);
        return
    end
end
fclose(fid);

% truncate variables
for isens= 1:nsens
    dat.(dat.sensorlist.name{isens})(icycle+1:end) = [];
end