### Code adapted from https://github.com/OpenMendel/SnpArrays.jl ###

"""
    makestream(filepath)

Open a file with an appropriate codec determined by file suffix.
"""
function makestream(filepath, args...; kwargs...)
    io = open(filepath, args...; kwargs...)
    if endswith(filepath, ".gz")
        codec = isreadonly(io) ? GzipDecompressor() : GzipCompressor()
    elseif endswith(filepath, ".zlib")
        codec = isreadonly(io) ? ZlibDecompressor() : ZlibCompressor()
    elseif endswith(filepath, ".zz")
        codec = isreadonly(io) ? DeflateDecompressor() : DeflateCompressor()
    elseif endswith(filepath, ".xz")
        codec = isreadonly(io) ? XzDecompressor() : XzCompressor()
    elseif endswith(filepath, ".zst")
        codec = isreadonly(io) ? ZstdDecompressor() : ZstdCompressor()
    elseif endswith(filepath, ".bz2")
        codec = isreadonly(io) ? Bzip2Decompressor() : Bzip2Compressor()
    else
        return io
    end
    TranscodingStream(codec, io)
end

function makestream(f::Function, args...)
    io = makestream(args...)
    try
        f(io)
    finally
        close(io)
    end
end

function parse_datapath(datapath)
    datadirname = splitdir(datapath)[2]
    temp = split(datadirname, '_')[3:end]
    parse.(Int,[temp[2], temp[4], temp[6], temp[8]])
end

function parse_genodir_name(genodir)
    genodirname = splitdir(genodir)[2]
    temp = split(genodirname, '_')[3:end]
    temp1 = parse.(Int,[temp[i] for i in 2:2:length(temp) - 1])
    rho = parse(Float64, temp[end])
    temp1, rho
end

function get_N(datapath::AbstractString)
    @assert isfile("$datapath/G1.fam") || isfile("$datapath/G.fam") "Genotype file name is not valid."
    srcfamfile = isfile("$datapath/G1.fam") ? "$datapath/G1.fam" : "$datapath/G.fam"
    srcN = makestream(srcfamfile) do stream
        countlines(stream)
    end
    srcN
end

function get_Mₖ(k::Integer, datapath::AbstractString)
    if k == 1
        srcbimfile = isfile("$datapath/G1.bim") ? "$datapath/G1.bim" : "$datapath/G.bim"
    else
        srcbimfile = "$datapath/G$k.bim"
    end
    srcMₖ = makestream(srcbimfile) do stream
        countlines(stream)
    end
    srcMₖ
end

get_K(genodir::AbstractString)=filter(contains(r"G\d+.bed"), readdir(genodir)) |> length

function get_avg_maf(k::Integer, datapath::AbstractString)
    src_maf_file = isfile("$datapath/MAF$k.txt") ? "$datapath/MAF$k.txt" : "$datapath/MAFs$k.txt"
    maf = makestream(src_maf_file) do stream
        maf = 0.0
        n = 0
        for line in eachline(stream)
            maf += parse(Float64, split(line)[1])
            n += 1
        end
        maf / n
    end
    maf
end
