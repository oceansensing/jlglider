using Glob, DataFrames, CSV, Dates, Missings
using TranscodingStreams, CodecZlib
import GZip

project = "NORSE"
deploydate = "20220818"
suffix = "test"
glidername = "sea064"
mission = "34"

yoid = 1

dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "nav/logs/";

missionroot = glidername * "." * mission;
gliroot = missionroot * "." * "gli.sub.";

datapathin = navdir * gliroot * string(yoid) * ".gz"; 
datapathout = navdir * gliroot * string(yoid) * ".csv";

#decompressed = transcode(GzipDecompressor, open(datapath));
#input = open(datapathin, "r");
#output = open(datapathout, "w");
#codec = GzipDecompressor();
#stream = TranscodingStream(codec, open(datapathin, "r"))

#stream = GzipDecompressorStream(input, stop_on_end=true);
#a = read(stream)

fh = GZip.open(datapathin)
a = readlines(fh)
close(fh)

