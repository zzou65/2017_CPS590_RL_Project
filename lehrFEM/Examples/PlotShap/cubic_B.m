% Run script for generating cubic bubble
  
%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    Cubic_Bubble([0 0;5 0;0 1]);

    print('-depsc', 'cubic.eps');
    !gv cubic.eps &