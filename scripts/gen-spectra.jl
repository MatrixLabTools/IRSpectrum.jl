#!/usr/bin/env julia

# Calculate spectra and save it to CSV file

using IRSpectrum
using CSV
using DataFrames
using ArgParse

function parse_commandline()
    s = ArgParseSettings(
        description = "Calculate spectra from CP2K dipole trajectory files and save it to CSV file."
    )

    @add_arg_table! s begin
        "--max_freq", "-f"
            help = "maximum frequency in wavenumbers"
            arg_type = Float64
            default = 4000.0
        "--time_step", "-t"
            help = "time step in femtoseconds"
            arg_type = Float64
            default = 0.5
        "--save_file", "-s"
            help = "output file"
            required = true
        "--use_filename"
            help = "save based on file name"
            action = :store_true
        "files"
            help = "trajectory files"
            nargs = '+'
            required = true
    end

    return parse_args(s; as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args: ")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    sp = DataFrame()
    tstep = parsed_args[:time_step]
    max_freq = parsed_args[:max_freq]
    sfile = parsed_args[:save_file]
    use_file_names = parsed_args[:use_filename]
    for (i,f) in enumerate(parsed_args[:files])
        s = spectrum(f; tstep=tstep, maxfreq=max_freq)
        if i == 1
            sp[!, :wavenumber] = s["wavenumber"]
        end
        if use_file_names
            name = split( splitpath(f)[end], ".")[begin]
            sp[!, Symbol(name)] = s["absorption"]
        else
            sp[!, Symbol("traj-$i")] = s["absorption"]
        end
    end
    CSV.write(sfile, sp)
end

main()