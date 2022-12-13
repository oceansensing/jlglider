function sensorlist = readcacheslocum(fname)

% Read Slocum glider sensor list from cachefile. Output as structure
% suitable for use in readslocumbd.m.
%
% Robert Todd, 29 April 2014

% parameters
maxsens = 3000;

% initialize
sensorlist.transmitted = false(maxsens,1);
sensorlist.sensnum = uint16(zeros(maxsens,1)); % will be 0:nsens-1
sensorlist.indexnum = int16(zeros(maxsens,1)); % could be -1 if not transmitted
sensorlist.nbytes = uint8(zeros(maxsens,1)); % 1 = 1 byte integer, 2 = 2 byte integer, 4 = 4 byte float, 8 = 8 byte double
sensorlist.name = cell(maxsens,1);
sensorlist.unit = cell(maxsens,1);

% read file
fid = fopen(fname);
nsens = 0;
while 1 % break loop at end of file
    line = fgetl(fid);
    if ~ischar(line), break, end
    nsens = nsens+1;
    x = textscan(line,'%*s %s %f %f %f %s %s');
    if strcmp(x{1}{1},'T')
        sensorlist.transmitted(nsens) = true;
    end
    sensorlist.sensnum(nsens) = uint16(x{2});
    sensorlist.indexnum(nsens) = int16(x{3});
    sensorlist.nbytes(nsens) = uint8(x{4});
    sensorlist.name{nsens} = x{5}{1};
    sensorlist.unit{nsens} = x{6}{1};
end
fclose(fid);

% truncate sensorlist
sensorlist.transmitted(nsens+1:end) = [];
sensorlist.sensnum(nsens+1:end) = [];
sensorlist.indexnum(nsens+1:end) = [];
sensorlist.nbytes(nsens+1:end) = []; 
sensorlist.name(nsens+1:end) = [];
sensorlist.unit(nsens+1:end) = [];