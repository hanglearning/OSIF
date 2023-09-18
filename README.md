This version of Open Source Impedance Fitter (OSIF) is modified by Hang Chen, Maria Van Venrooy, Lalith Teja Nagidi, Nihaal Chowdary and Emma Cao during the University of Delaware 2023 DSI Hackathon for the Chemours company in Newark, Delaware that supports the following serialization functionalities compared to the orignla OSIF software (https://github.com/NREL/OSIF)

1. Indicate lower bound, upper bound, and value internval of model parameters of the selected models
2. Automated evaluation and selection of best fitting parameters and the model based on cost value fitting to the base datafile
3. Serialize the parameters and the model having the lowest cost to the rest of the data files

Our project extends the National Renewable Energy Laboratory's Open-Source Impedance Fitter (OSIF), enabling it to autonomously optimize model parameters for electrochemical impedance spectra from Proton-exchange membrane fuel cells in hydrogen/nitrogen environments. These parameters align with quasi-transmission line, one-dimensional linear diffusion, or spherical diffusion models. The result is an upgraded software version featuring a user-friendly interface (UI) for streamlined functionality. It excels at automatically applying optimal model parameters to additional datasets, significantly reducing the need for manual trial and error. This innovation stands to save valuable time for chemical engineers and researchers. Our software empowers swift modeling and analysis of Proton-exchange membrane fuel cell impedance data, enhancing research and development in electrochemistry. Its accessibility through a user-friendly interface broadens its utility, facilitating broader exploration and adoption of these critical technologies. Our project represents a crucial advancement in fuel cell technology research and engineering efforts.

Here's the introduction of the original OSIF software copied from https://github.com/NREL/OSIF

Open Source Impedance Fitter (OSIF) is a program that allows the user to fit electrochemical impedance spectra of Proton-exchange membrane fuel cells collected under a hydrogen/nitrogen (anode/cathode) environment to an accepted quasi-transmission line model (1), or one dimentional linear diffusion or spherical diffusion models (2).

To use OSIF, the following non-default python packages will have to be downloaded:
  - xlrd
  - matplotlib
  - numpy
  - scipy
  
The default python version the program is written for is 3.x, though it can be used in python 2.x (see in "How to Use OSIF" file in the repository)

The user must have impedance data in one of the following formats (see example data folder in repository):
  - An excel spread sheet with columns as, Index, Frequency (Hz), Z' (Ω), -Z'' (Ω), Z (Ω), -Phase (°), Time (s)
  - A tab delimited text file with a header depicting the columns as Index, Frequency (Hz), Z' (Ω), -Z'' (Ω), Z (Ω), -Phase (°), Time (s) 
  
Where in Z' is the real part of the measured impedance, Z'' is the imaginary part of the measured impedance, and Z is the modulus of the measured impedance.
  
  
For instructions on how to, install python on a windows or mac platform (in order to use OSIF) and running OSIF as well as more information on how OSIF works see the "How to Use OSIF" file in the repository.



(1) Setzler, B. P., & Fuller, T. F. (2015). A Physics-Based Impedance Model of Proton Exchange Membrane Fuel Cells Exhibiting Low-Frequency Inductive Loops. Journal of The Electrochemical Society, 162(6), F519-F530. doi:10.1149/2.0361506jes
(2) Jean-Paul Diard (BioLogic), Bernard Le Gorrec (UJFG-INPG), C Montella. Handbook of Electrochemical Impedance Spectroscopy. DIFFUSION IMPEDANCES.
https://www.researchgate.net/publication/342833389_Handbook_of_Electrochemical_Impedance_Spectroscopy_DIFFUSION_IMPEDANCES
