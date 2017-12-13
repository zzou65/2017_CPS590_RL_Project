% Run script for generating pyramid by the comination of LFE shap function
  
%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Initialize constant
    
    NPts = 10;
    
    % Plot out the pyramid
    
    plot_Shap_LFE(NPts,1,[0 0;1 0;0 1])
    plot_Shap_LFE(NPts,1,[0 0;0 -1;1 0])
    plot_Shap_LFE(NPts,1,[0 0;-1 0;0 -1])
    plot_Shap_LFE(NPts,1,[0 0;0 1;-1 0])
    
    set(gca,'CameraPosition', [1 -.5 1]*50)
    
    % Output .eps files
    
    FileName = ['Shap_Pyramid3D_LFE.eps'];
    print('-depsc', FileName);
    system(['gv ' FileName ' &']);
    
    