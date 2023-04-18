# example loading of .mat file produced by odas_p2mat.m
# 2023-04-18 gong@vims.edu

using Glob, MAT

datadirJM = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/mr1000g/"
datadirLBE = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/mr1000g/"

datadir = datadirJM;
datafiles = Glob.glob("*.mat", datadir);

file = matopen(datafiles[1])
a = read(file, "Gnd")
close(file)

