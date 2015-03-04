HHG_Phase-matching
=================

Takes the equation from Constant et al. PRL 82, 1668 - 1671 (1999) and implements them.

## Brief file description
### /src
`atom.py` contains useful atomic parameters to calculate things like dispersion at the driving wavelength and harmonics. It also contains parameters to perform ADK ionization rate calculations.

`ADK.py` contains two version of an ADK ionization rate calculation. It can provide the ionizationrate at a given intensity or integrated rate over a pulse.

`laser.py` contains laser parameters and two different types of pulses. Sin ** 2 pulse and gaussian pulse. 

`gas.py` contains gas jet parameters.

`sf.py` and `/sf` contain atomic scattering factors provided by the LBL. The values for f1 below 30 eV were calculated using the Kramers-Kronig relations and the measured values of f2.

`phasematching.py` is where most of the work is done. It can provide the harmonic yield and parameters like the coherence length as a function of time across a pulse or an integrated yield over the pulse. It takes all the relevent simulation parameters.

### /test
Most of the tests runs produce matplotlib plots of the results. 
`
`testrun.py`, `testrun2.py`, and `testrun3.py` contain example runs and compares them to other results.

`i_scan.py` is an intensity scan of the laser.

`p_scan.py` is a pressure scan.

`w_scan.py` is a laser focusing spot size scan.

`pulse_scan.py` is a pulse duration scan. 

### /images
Contains various outputs of the code to compare with known experimental results. The file name gives a hint as to what, but probably only means anything to me, sorry...



