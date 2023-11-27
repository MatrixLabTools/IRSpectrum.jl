module IRSpectrum

using FFTW

export autocorrelation,
       read_cp2k_dipoles,
       spectrum


"""
read_cp2k_dipoles(fname) -> Matrix{Float64}

Read dipole trajectory file created by CP2K and return dipolemoments as
a two dimensional array.
"""
function read_cp2k_dipoles(fname::AbstractString)
    out = Float64[]
    open(fname,"r") do file
        re = r"^\s*X=\s+(?<X>-?\d+.?\d+[Ee]?[+\-]?\d+).*\s+Y=\s+(?<Y>-?\d+.?\d+[Ee]?[+\-]?\d+).*\s+Z=\s+(?<Z>-?\d+.?\d+[Ee]?[+\-]?\d+)"
        for line in eachline(file)
            m = match(re, line)
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
autocorrelation(μ::AbstractMatrix) -> Vector

Calculates autocorrelation function using [Wiener-Khinchin Theorem](http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html).
"""
function autocorrelation(μ::AbstractMatrix)
    tmp = [ fft(abs2.(ifft(r))) for r in eachrow(μ)]
    return real.(reduce(+,tmp))
end


"""
spectrum(μ::AbstractMatrix; tstep=0.5, maxfreq=4000) -> Dict

Calculates spectrum from dipolemoment trajectory. Calculates autocorrelation function and
its Fourier transform.

# Arguments
- `μ::AbstractMatrix` :  data where autocorrelation funtion is calculated for each row

# Keywords
- `tstep` : timestep for trajectory in fempto seconds
- `maxfreq` : maximum frequency in wavenumbers for resulting spectrum

# Returns
- `Dict` : with fields `"absorption"` and `"wavenumber"`.
"""
function spectrum(μ::AbstractMatrix; tstep=0.5, maxfreq=4000)
    l = size(μ)[2]
    lh = Int(floor(l/2))
    acor = autocorrelation(μ)
    s = abs.(fft(acor))[1:lh]
    # freq = n/NΔ and wavenumber = freq/100c
    w = [ n/(l*tstep*1E-15*299792458E2) for n in 1:lh ]
    i = w .<= maxfreq
    return Dict("absorption"=>s[i], "wavenumber"=>w[i])
end


"""
spectrum(fname::AbstractString; tstep=0.5, maxfreq=4000) -> Dict

Read dipolemoments from CP2K dipole file `fname` and calculate spectrum.
"""
function spectrum(fname::AbstractString; tstep=0.5, maxfreq=4000)
    μ = read_cp2k_dipoles(fname)
    return spectrum(μ; tstep=tstep, maxfreq=maxfreq)
end

end # module
