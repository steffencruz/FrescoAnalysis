# FrescoAnalysis
DWBA Analysis Companion

- Simple to use Fresco Manager
- Uses FRESCO template files and neutron state information to generate a DWBA calculation and fit it to data to extract experimental spectroscopic factors
- Can also be used to generate elastic scattering files to normalize data to DWBA

  
  
# Search Mode

- Calculates SF for all possible neutron states and returns lowest chi2 fit

# Scan Mode

- Calculates SF & Sig for all possible neutron states and produces a graphs of these values as a function of ExcHi and fits this to a cumulative function which indicates the best SF.

# Description of variables

- REAC - "dp", "pp", "dd" selects the reaction

** neutron state information **  
- EXC - Excitation energy [MeV]
- BE  - Binding energy of neutron [MeV]
- JF  - Total angular momentum of final state, coupling [s1/2][nlj]
- J0  - Total angular momentum of orbital that neutron is placed in
- L0  - Orbital angular momentum of orbital that neutron is placed in
- N0  - Principal quantum number of orbital that neutron is placed in
- SF  - Spectroscopic Factor
- CS - Total cross section [with SF appplied] 

** required to determine spectroscopic factor error **
 - NORM - Normalization Used
 - GAM  - Gamma gate used 
 
** these are taken from the Results_*Exc*.root file **
- ExcHi - Upper excitation energy window, see ScanEnergyWindow

# Verbose : -
 - 0   =   No Printing 
 - 1   =   Only High level summary eg. CheckFeeding
 - 2   =   Standard printing
 - 3   =   Print everything


# Notes

Further work : -
- Add (d,t) reaction
- Copy files and make directories?
- Speed up the process of doing calculations. Scripts.
   
Known bugs :-
