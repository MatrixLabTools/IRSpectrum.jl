# IRSpectrum.jl

This package is meant to calculate Vibration spectrum from molecular dynamics
data. Mainly by using dipolemoment trajectory to calculate IR-spectrum or
velocity trajectory to calculate power spectrum.

Main user interface is `spectrum`-function that will calculate spectrum from
given data that can be either dipole moments or velocites. Autocorrelation
function is calculated with
[Wiener-Khinchin Theorem](http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html)
and it is thus scaling with NlogN.

There are also special readers for CP2K dipole output. You can either read
dipoles to a matrix and then calculate spectrum or just call `spectrum` with
file name to get IR-spectrum right away.

### To calculate IR-spectrum use

```
using IRSpectrum

# tstep is time step in fs
# maxfreq is maximum frequency in wavenumbers
spectrum(path_to_cp2k_dipole_traj_file; tstep=0.5, maxfreq=4000)
```
