# Modeling AIR-1, ECT-2 and Myosin

This folder contains Matlab codes for the publication
"A minimal mathematical model for polarity establishment and centralsplindlin-independent cytokinesis," 
by O. Maxian, K. Longhini, and M. Glotzer, Aug. 2024.

There are three files:
* DiffusionEmbryo.m - this is the file that computes the AIR-1 signal on the cell (ellipse) boundaries from the centrosome positions (Appendix A of paper). 
* MyosinEct2.m - this file runs the dynamics of myosin and Ect-2 forward in time using the equations and parameters in Appendix B of the paper. 
Running it will generate Fig. 4(c). 
* MyosinEct2_Stability.m - performs the stability analysis on the model in Fig. S4. 

In addition, the folder "Data" contains the saved simulated AIR-1 signals as follows:
* LonghiniData.mat has a 9 x 4 array of centrosome positions and AIR-1 signals from experimental measurements in cytokinesis (Longhini and Glotzer, 2022).
* AIR1Diffusion.mat contains the simulated AIR-1 signals during cytokinesis across each of the nine embryo treatments as a function of arclength. 
* PolarizationAllDistances.mat and PolarizationDistances.mat contains the simulated AIR-1 signals during polarization for different centrosome distances. 
* PolarizationKoff.mat contains the simulated AIR-1 signals for different values of koff (see Fig. S1)
