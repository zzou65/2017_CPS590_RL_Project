% Run script for uniform similarity 

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Initialize constant
    
    MM = 50;
    tol = 1e-5;
    H0 = 0.01;
    BBOX = [0 0;1 1];
    HHANDLE = @h_uniform;    
    DHANDLE = inline(['dist_isect(dist_isect(dist_circ(x,[0 0],1),dist_circ(x,[1 0],1)),' ...
                                 'dist_rect(x,[0 0],1,1))'],'x');    
    FIXEDPOS = [0 0;1 0;.5 3^.5/2];  
    DISP = 0;
    
    % Preallocate memory
    
    % Generate the mesh
    
    Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);

    % Compute the function
    
    nCoordinates = size(Mesh.Coordinates,1);
    U = zeros(nCoordinates,1);
    
    for i=1:nCoordinates
        
        vert = Mesh.Coordinates(i,:);
        px = vert(1);
        py = vert(2);
        a = norm(vert);
        b = norm(vert-[1 0]);
        if (norm(py) < tol )
            U(i) = log10(MM);
        else
            U(i) = log10(min((a+b+1)/py,MM));
        end
        
    end
    
    % Plot the solution
    
    plot_LFE(U,Mesh);
    colorbar;
    title('\bff(x,y)=log_{10}((radius of inscribed circle)^{-1})')
    xlabel('\bfp_{x}')
    ylabel('\bfp_{y}')
    % Output .eps files
    
    print('-depsc', 'uni_insc.eps');
    !gv uni_insc.eps &
    
    % Clear memory
    
    clear all;