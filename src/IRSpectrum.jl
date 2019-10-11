module IRSpectrum

using FFTW

export read_cp2k_dipoles, autocorrelation, spectrum


"""
read_cp2k_dipoles(fname)

Read dipole trajectory file created by CP2K and return dipolemoments as
a two dimensional array.
"""
function read_cp2k_dipoles(fname)
    out = Float64[]
    open(fname,"r") do file
        for line in eachline(file)
            m = match(r"X=\ +(?<X>-?\d+.?\d+)\ +Y=\ +(?<Y>-?\d+.?\d+)\ +Z=\ +(?<Z>-?\d+.?\d+)", line)
            if m !== nothing
                push!(out,parse(Float64,m[:X]))
                push!(out,parse(Float64,m[:Y]))
                push!(out,parse(Float64,m[:Z]))
            end
        end
    end
    if length(out)%3 == 0
        return reshape(out, (3, Int(length(out)/3)))
    else
        error("Error in reading file")
    end
end


"""
autocorrelation(μ::AbstractMatrix)

Calculates autocorrelation function using [Wiener-Khinchin Theorem](http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html).
"""
function autocorrelation(μ::AbstractMatrix)
    tmp = [ fft(abs2.(ifft(r))) for r in eachrow(μ)]
    return real.(reduce(+,tmp))
end


"""
spectrum(μ::AbstractMatrix; tstep=0.5, maxfreq=4000)

Calculates spectrum from dipolemoment trajectory.
`tstep` is timestep for trajectory in fempto seconds and
`maxfreq` is maximum frequency in wavenumbers for resulting spectrum.

Return `Dict` with field `"absorption"` and `"wavenumber"`.
"""
function spectrum(μ::AbstractMatrix; tstep=0.5, maxfreq=4000)
    l = size(μ)[2]
    lh = Int(floor(l/2))
    acor = autocorrelation(μ)
    s = abs.(fft(acor))[1:lh]
    w = [ n/(l*0.5E-15*299792458E2) for n in 1:lh ]
    i = w .<= maxfreq
    return Dict("absorption"=>s[i], "wavenumber"=>w[i])
end


"""
spectrum(fname::AbstractString; tstep=0.5, maxfreq=4000)

Read dipolemoments from file `fname` and calculate spectrum.
"""
function spectrum(fname::AbstractString; tstep=0.5, maxfreq=4000)
    μ = read_cp2k_dipoles(fname)
    return spectrum(μ, tstep=tstep, maxfreq=maxfreq)
end

end # module
