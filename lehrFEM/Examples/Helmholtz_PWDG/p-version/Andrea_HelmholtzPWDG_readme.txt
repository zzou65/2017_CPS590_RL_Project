12.5-10.12.2009 Andrea
Files created / modified in LehrFEM-PWDG for experiments for paper p-convergence, April/May 2009, changed a bit in December

LehrFEM/Examples/Helmholtz_PWDG/P-Conv:
_conv_p_and.m		main file created to plot the p-convergence of the errors
_conP_omegas.m		file to plot the convergence with different omega, on the same mesh or scaling meshes
_plot_cornerwave.m	plot the exact solution

_LInfErr_PWDG.m			compute error in L^2 and also L^\infty
_MultiNormErr_PWDG.m		compute separated errors in L^2, broken H1 s.n., L^2(jumps)
(this two could be put in   LehrFEM/Lib/Errors/ )


June 2010 (for Jay)
_createPWDGdata.m	solve the equation with different parameters and create a .mat file with the errors
_conv_from_data.m       various plot with the data from PWDG_error_data.mat
_PWDG_error_data.mat	data created by createPWDGdata.m


_figs\			subfolder with figures put in the paper (except the two from Ralf)
	edgesing\	singularity on an edge instead of corner, only minor differences
	ord\		plot (boh?) with lines for slopes for Ilaria
