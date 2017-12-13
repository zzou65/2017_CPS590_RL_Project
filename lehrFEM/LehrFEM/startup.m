% Startup script for LehrFEM library
  
% Copyright 2005-2005 Patrick Meury, Kah-Ling Sia & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Print header
  
  fprintf('\n');
  fprintf('   ###########################################\n');
  fprintf('   ###                                     ###\n');
  fprintf('   ### LehrFEM - 2D finite element toolbox ###\n');
  fprintf('   ###                                     ###\n');
  fprintf('   ###########################################\n');
  fprintf('\n');
  
  % Add assembly routines
  
  fprintf('     Add assembly routines                     ');
  try
    addpath([pwd '/Lib/Assembly']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end

  % Add element compuation routines
  
  fprintf('     Add element computation routines          ');
  try
    addpath([pwd '/Lib/Element']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add mesh generation routines routines
  
  fprintf('     Add mesh generation routines              ');
  try
    addpath([pwd '/Lib/MeshGen']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');
  end
  
  % Add plotting routines
  
  fprintf('     Add plotting routines                     ');
  try
    addpath([pwd '/Lib/Plots']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add quadrature rules
  
  fprintf('     Add quadrature rules                      ');
  try
    addpath([pwd '/Lib/QuadRules']); 
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add error estimators
  
  fprintf('     Add error estimators                      ');
  try
    addpath([pwd '/Lib/ErrEst']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add discretization errror routines
  
  fprintf('     Add discretization errors                 ');
  try
    addpath([pwd '/Lib/Errors']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add finite element functionals
  
  fprintf('     Add functionals                           ');
  try
    addpath([pwd '/Lib/Functionals']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add solver routines
  
  fprintf('     Add solver routines                       ');
  try
    addpath([pwd '/Lib/Solvers']);
    addpath([pwd '/Lib/Solvers/AMGToolbox/matlab_amg_release_Oct2002']);
    addpath([pwd '/Lib/Solvers/AMGToolbox/matlab_util_release_Oct2002']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Time stepping routines
  
  fprintf('     Add time stepping routines                ');
  try
    addpath([pwd '/Lib/TimeStep']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add tool finite element tools 
  
  fprintf('     Add tools                                 ');
  try
    addpath([pwd '/Lib/Tools']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');
  end
  
  % Add error distribution routines
  
  fprintf('     Add error distribution routines           ');
  try
    addpath([pwd '/Lib/ErrorDistr']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add interpolation routines
  
  fprintf('     Add interpolation routines                ');
  try
    addpath([pwd '/Lib/Interp']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');
  end
  
  % Add limiter routines
  
  fprintf('     Add limiter routines                      ');
  try
    addpath([pwd '/Lib/Limiters']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  % Add m2htlm-Tool routines
  
  fprintf('     Add m2htlm-Tool routines                  ');
  try
    addpath([pwd '/Lib/Tools/m2html']);
    fprintf('[  OK  ]\n');
  catch
    fprintf('[FAILED]\n');  
  end
  
  fprintf('\n')