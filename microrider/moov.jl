
# the moov function moves files matching a specified string pattern from fromdir to todir
# Donglai Gong, 2024-01-14 : initial write
#
# example usage:
#stringpattern = "*.txt";
#fromdir = "director1";
#todir = "directory2";
#moov(stringpattern, fromdir, todir)

using Glob
function moov(stringpattern::String, fromdir::String, todir::String)
    frompath = readdir(Glob.GlobMatch(stringpattern), fromdir);
    topath = joinpath.(todir, basename.(frompath));
    mv.(frompath, topath, force=true);
end
