#include("run_seaexplorer");

x1 = jm.lon;
y1 = jm.lat;
z1 = jm.z;
c1 = jm.ctemp;

z1ind = findall( -30 .<= z1 .<= -10);
pind1 = z1ind;

x2 = lbe.lon;
y2 = lbe.lat;
z2 = lbe.z;
c2 = lbe.ctemp;

z2ind = findall( -30 .<= z2 .<= -10);
pind2 = z2ind;