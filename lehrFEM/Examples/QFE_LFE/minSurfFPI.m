% Run script for minimal surface problem with the use of Fixed Point
% Iteration

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

% Initialize constants

    TOL = 0.0000001;
    MAXIT = 500;
    NREFS = 3;
    PlotFreq = 7;
    % Data for circular mesh

        GD_HANDLE = @(x,varargin)x(:,1);
        U0 = @(x,varargin)x(:,1)+x(:,2).^2+100*norm(x(:,1))*norm(x(:,2))*(1-sqrt(x(:,1).^2+x(:,2).^2))*norm(sin(pi*x(:,1)./50));
        DHANDLE = @dist_circ;
        C=[0,0];
        R=1;
    % Data for squared mesh

        %GD_HANDLE = @(x,varargin)x(:,1);
        %U0 = @(x,varargin)x(:,1)+10*(x(:,1)-1).*x(:,1).*x(:,2).*(x(:,2)-1);

    % Data fpr LShaped Mesh

        %GD_HANDLE = @(x,varargin)2;                       % Boundary data
        %U0 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2)); % Initial guess of solution

    % Initialize mesh
    %Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
    Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
    %Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1:-1:-4;
    for i = 1:NREFS,
      %Mesh = refine_REG(Mesh);
      Mesh = refine_REG(Mesh,DHANDLE,C,R);
    end
    Mesh = add_Edge2Elem(Mesh);

    % Interpolate initial data

    M = assemMat_LFE(Mesh,@MASS_LFE);
    L = assemLoad_LFE(Mesh,P3O3(),U0);
    U = M\L;
    U0 = U;
    plot_LFE(U,Mesh);
    colorbar;
    title('{\bf Initial guess}');
    xlabel(['{\bf ' int2str(0) ' iterations}']);
    drawnow();    

    err = 2 * TOL;
    it = 0;per = 0;
    while err > TOL,
       it = it + 1;
       if it > MAXIT,
           error('Too many iterations');
       end
       U_old = U; % Keep this for termination criteria

       A = assemMAT_minSurf(U,Mesh,@STIMA_minSurf);
       L = zeros(size(Mesh.Coordinates,1),1);

       % Incorporate Dirichlet boundary data

       [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
       L = L - A*U;

       % Solve the linear system

       U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);

       % Check if we are done

       err = norm(U(FreeDofs)-U_old(FreeDofs))/norm(U0(FreeDofs));

       if rem(it,PlotFreq) == 0,
           plot_LFE(U,Mesh);
           colorbar;
           title('{\bf Approximate solution}');
           xlabel(['{\bf ' int2str(it) ' iterations}']);
           drawnow();

       end

    end

% Final plot (fancy plot)

plot_LFEfancy(U,Mesh);