%Algebraic Multigrid Package
%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% These programs were prepared by the Regents of the University of
%% California at Los Alamos National Laboratory (the University) under
%% contract W-7405-ENG-36 with the U.S. Department of Energy (DOE) and
%% under Contract KC-07-01-01 with the Department of Energy(DOE),
%% Office of Science, Mathematical, Information, and Computational Sciences,
%% Applied Mathematical Sciences Program.  All rights in these programs
%% are reserved by the DOE and the University.

%% Permission is granted to the public to copy and use this software
%% without charge, provided that this Notice and the statements of 
%% authorship are reproduced on all copies. Neither the U.S. government 
%% nor the University makes any warranty, express or implied, or assumes 
%% any liability or responsibility for the use of this software.

%% AMG code written by Menno Verbeek, in collaboration with Jane Cullum
%%  and Wayne Joubert.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Setup
%   AMGSetup		 - AMG setup phase
%   AMGDefaultOptions    - Return defautl options
%   AMGSelectCoarseGrid  - Selects the course grid points (Ruge-Stueben)
%   AMGMakeInterpolation - Make the interpolation matrix P
%   Ruge_Stueben   	 - AMG interpolation matrix using Ruge-Stueben alg.
%   GetDistances    	 - Get distance matrix from the grid (2D only)
%   readall       	 - Read matrix output from lamg and store in AMG data 
%   mesh2coord    	 - Convert mesh data to grid variable coordinates
%   propagategrid 	 - Propagate grid data down the levels
%   standardgrid   	 - Make grid info for standard 2-D 5-point stencil
%                          structure
%
%Operators
%   AMGVcycle     	 - Preform an AMG V-cycle
%   AMGBlindErrorOperator- Error operator for the AMG V-cycle SKS=(1-Vcycle*A) 
%   AMGBlindS       	 - Error operator for the AMG smoother S=(1-Mpre*A)
%   AMGBlindK       	 - Error operator for the AMG coarse gird correction K
%
%Diagnostics 
%   plotprepost    	 - Do full eigen analisys of AMG problem
%   ComputeProperties    - Compute the properties of an AMG setup
%   AMGconv       	 - Get maximum eigenvalue of overall AMG erorr operator
%   changebasis    	 - Change AMG data structure to P and R basis
%
%Plotting tools
%   AMGPlotGird     	 - Plot the mesh, CF-split and matrix connectivity
%   AMGPlotVec      	 - Plot a vector on the grid
%   PlotTable     	 - Plot table using pcolor plot
%   plotgrids       	 - Plot AMG grid points for the different levels
%   plotmesh      	 - Plot mesh data
%   showlocal      	 - Attempt to plot local matrix connenctions in graph
%
%Test Runs
%   RunTests        	 - Run tests with different matrices and options
%   FinishTests   	 - Finish RunTests run or fill empty spots in RunTests 
%                          data structure
%   MergeResults  	 - Merge RunTests results
%   ReorderResults  	 - Reoder Options in AMG RunTests results
%   ResultSortNames	 - Reoder Names in AMG RunTests results
%   SaveResults		 - Save results struct to file
%   LoopOptions     	 - Generate options struct array for RunTests
%   LoopOptionsRS     	 - Generate Ruge-Stueben options
%   LoopOptionsRSFill    - Generate Ruge-Stueben options with fill
%   LoopOptionsEig   	 - Generate eigenvector options
%   LoopOptionsBig	 - Generate limited set of options for big matrices
%   GetProp		 - Extract property matrix from AMG RunTests Prop struct
%   OptionsTable   	 - Make a .tex and .ps table of an options struct array
%   Table2Tex      	 - Make latex and .ps table output
%
%NOTE: The AMG code needs helper routines form the AMGutil directory.
%


%Helper routines
%   gsai          	 - Compute generalised sparse approximate inverse
%   ainv_explicit 	 - Get explicit approximate inverse using ainv
%   ainv          	 - Compute approximate factorization of inverse of 
%   arnoldi       	 - Block Arnoldi recursion
%   deflate         	 - Part of block Arnoldi routine "arnoldi"
%   iorder          	 - Part of block Arnoldi routine "arnoldi"
%   myarnoldi      	 - Very simple Arnoldi algorithem
%                          SPD matrix
%
