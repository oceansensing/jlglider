using Glob, DataFrames, CSV, Dates, Missings
import TranscodingStreams, CodecZlib
import GZip

datadir = dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";

transcode(GzipDecompressor, open()