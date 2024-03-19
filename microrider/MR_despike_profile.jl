using FFTW, DSP

figoutdir = "/Users/gong/GitHub/jlglider/microrider/figures/";

function mr_despike(sh1_hp, spikethreshold)
    # Applying the lowpass filter
    sh1_lp = DSP.filtfilt(digital_filter_lp, sh1_hp);

    # find the indices in the cast that exceed the spike threshold (typically 8)
    badind = findall(sh1_hp ./ sh1_lp .> spikethreshold);

    # find the indices in the cast that stays below the spike threshold (typically 8)
    #godind = findall(sh1_hp ./ sh1_lp .<= spikethreshold);

    # extract only those indices that are within the depth range of interest (i.e. avoid surface and bottom inflections)
    bzind = badind;
    #gzind = godind;
    #bzind = intersect(badind, zinddn);
    #gzind = intersect(godind, zinddn);

    # Apply the highpass filter
    sh1_hp = abs.(DSP.filtfilt(digital_filter_hp, sh1_hp));

    N = 20;
    for ii = 1:length(bzind)
        #display(ii)
        bseg = Int64.(bzind[ii]-N/2 : bzind[ii]+N); # bad segment due to spike
        rind = [collect(bseg[1]-256 : bseg[1]-1) ; collect(bseg[end]+1 : bseg[end]+256)]; # replacement value indices
        rind = rind[(rind .> 0) .&& (rind .<= length(sh1_hp))];
        bseg = bseg[(bseg .> 0) .&& (bseg .<= length(sh1_hp))];

        sh1_repval = mean(sh1_hp[rind]); # shear 1 replacement value for this spike segment
        sh1_hp[bseg] .= sh1_repval;
    end

    #sh1_hp = sh1_hp ./ maximum(sh1_hp);

    nbadind = length(findall(sh1_hp ./ sh1_lp .> spikethreshold)); 
    return sh1_hp, nbadind;
end

ii = 1
zedge = 10;

mrp = norse23mr[ii].mr;
mrpz = norse23mr[ii].z[:];

sh1_raw = mrp.sh1;
sh2_raw = mrp.sh2;

minzind = findall(mrpz .== minimum(mrpz))[1]
tbott = round(mrp.t_fast[minzind], RoundDown);
tinddn = findall(tbott .>= mrp.t_fast[:] .>= 0);
zindOuter = findall(mrpz[minzind]+zedge-2 .<= mrpz[tinddn] .<= -zedge+2);
zindInner = findall(mrpz[minzind]+zedge .<= mrpz[tinddn] .<= -zedge);

# Generate a sample signal (replace this with your actual signal)
fs = 512  # Sample rate (Hz)
#t = 0:1/fs:tbott  # Time vector
nyquist_rate = fs / 2
spikethreshold = 8.0;

# Design a high-pass filter
cutoff_frequency_hp = 0.1  # Cutoff frequency (Hz)
normalized_cutoff_hp = cutoff_frequency_hp / nyquist_rate
digital_filter_hp = digitalfilter(DSP.Highpass(normalized_cutoff_hp), DSP.Butterworth(1))

# Design a low-pass filter
cutoff_frequency_lp = 0.3; # Cuttoff frequency (Hz)
normalized_cutoff_lp = cutoff_frequency_lp / nyquist_rate
digital_filter_lp = digitalfilter(DSP.Lowpass(normalized_cutoff_lp), DSP.Butterworth(1))

global sh1_hp = abs.(DSP.filtfilt(digital_filter_hp, mrp.sh1[tinddn[zindOuter]]));
global sh2_hp = abs.(DSP.filtfilt(digital_filter_hp, mrp.sh2[tinddn[zindOuter]]));
#global sh1_hp = sh1_hp ./ maximum(sh1_hp);

global itr = 0;
nbadind = 1;
#while ((nbadind == 0) & (itr <= 2)) | (nbadind > 0)
for ii = 1:10
    global itr = itr + 1;
    global sh1_hp, nbadind = mr_despike(sh1_hp, spikethreshold);
    global sh2_hp, nbadind = mr_despike(sh2_hp, spikethreshold);
    #if itr <= 4
    #    nbadind = 1;
    #end
    display([itr nbadind])
end

zouter = mrpz[tinddn[zindOuter]];
touter = mrp.t_fast[tinddn[zindOuter]];
zind = findall(mrpz[minzind]+zedge .<= zouter .<= -zedge); 
zinner = zouter[zind];
tinner = touter[zind];

#Plots.plot(sh1_hp[zind])

GLMakie.closeall()
psize = (2000, 1000);
figsh = Figure(size = psize, fontsize = 32)

axsh1 = Axis(figsh[1, 1],
    title = "Shear 1 (HP, filtered)",
    xlabel = "Shear 1",
    ylabel = "Depth (m)",
)
GLMakie.scatterlines!(sh1_hp[zind], zinner, markersize = 3, linewidth = 0.5);
#GLMakie.xlims!(0, 0.001);
#GLMakie.ylims!(-12,-11.7);

axsh2 = Axis(figsh[1, 2],
    title = "Shear 2 (HP, filtered)",
    xlabel = "Shear 2",
    ylabel = "Depth (m)"
)
GLMakie.scatterlines!(sh2_hp[zind], zinner, markersize = 3, linewidth = 0.5);

axAx = Axis(figsh[1, 3],
    title = "Ax",
    xlabel = "Ax",
    ylabel = "Depth (m)"
)
GLMakie.scatterlines!(mrp.Ax[tinddn[zindOuter[zind]]], mrpz[tinddn[zindOuter[zind]]], markersize = 3, linewidth = 0.5);

axAx = Axis(figsh[1, 4],
    title = "Shear 1 (Raw)",
    xlabel = "Sh1",
    ylabel = "Depth (m)"
)
GLMakie.scatterlines!(sh1_raw[zind], zinner, markersize = 3, linewidth = 0.5);

GLMakie.Label(figsh[0, :], text = project * ": " * mission * " " * glider * " - " * basename.(mrp.fullPath[1:end-2]), fontsize = 50)

figsh

save(figoutdir * project * "_" * mission * "_" * glider * "_" * string(profileid[ii], pad = 4) * "_shear_profile_" * mrp.date * ".png", figsh)

#GLMakie.closeall()