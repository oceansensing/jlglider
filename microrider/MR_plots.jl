using GLMakie

figoutdir = "/Users/gong/GitHub/jlglider/microrider/figures/";

ip = 1;
mrp = norse23mr[ip].mr;
mrpz = norse23mr[ip].z;

zrange = [-150, 0];

#mrpz = gsw.gsw_z_from_p.(mrp.P_fast[:], 71.0, 0.0, 0.0);
ind = findall((zrange[2] .>= mrpz .>= zrange[1]) .& (mrp.t_fast[:] .< maximum(mrp.t_fast)/2));

xfast = mrp.t_fast[ind];
#yfast = mrp.P_fast[:];
yfast = mrpz[ind];
z1fast = mrp.sh1[ind];
z2fast = mrp.sh2[ind];

x = xfast;
y = yfast;
z = z1fast;

zmin = -0.02;
zmax = 0.02; 

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 3; # day intervals for plotting
ms = 6;
psize = (4000, 4000)


figP = Figure(; size = psize, fontsize = 64)
axP = Axis(figP[1, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " Shear 1",
    xlabel = "Time",
    ylabel = "Depth (m)"
)
GLMakie.plot!(x, y, color=z, colormap=:jet, markersize=40, colorrange=(zmin, zmax))
#ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(figP[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
figP
save(figoutdir * project * "_" * mission * "_" * glider * "_" * string(profileid[ip], pad = 4) * "_sh1.png", figP)
GLMakie.closeall()


figshear = Figure(size = psize, fontsize = 48)

axsh1 = Axis(figshear[1, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " sh1",
    xlabel = "Time",
    ylabel = "Shear 1",
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.sh1[ind], markersize = 3, linewidth = 0.5);

axsh2 = Axis(figshear[2, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " sh2",
    xlabel = "Time",
    ylabel = "Shear 2"
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.sh2[ind], markersize = 3, linewidth = 0.5);

axAx = Axis(figshear[3, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " Ax",
    xlabel = "Time",
    ylabel = "Ax"
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.Ax[ind], markersize = 3, linewidth = 0.5);

axAy = Axis(figshear[4, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " Ay",
    xlabel = "Time",
    ylabel = "Ay"
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.Ay[ind], markersize = 3, linewidth = 0.5);

axGradT1 = Axis(figshear[5, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " gradT1",
    xlabel = "Time",
    ylabel = "gradT1"
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.gradT1[ind], markersize = 3, linewidth = 0.5);

axGradT2 = Axis(figshear[6, 1],
    title = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]) * " gradT2",
    xlabel = "Time",
    ylabel = "gradT2"
)
GLMakie.scatterlines!(mrp.t_fast[ind], mrp.gradT2[ind], markersize = 3, linewidth = 0.5);

figshear
save(figoutdir * project * "_" * mission * "_" * glider * "_" * string(profileid[ip], pad = 4) * "_shear_ts.png", figshear)

#GLMakie.closeall()