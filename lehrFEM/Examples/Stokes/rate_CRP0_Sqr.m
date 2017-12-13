% Run script for the convergence rate for CR stokes solver on the unit 
% square domain.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear R_Mesh;
    
    % Initialize constant
    
    NREFS = 5;
    NU = 1;  
    F_Handle = @(x,varargin)[ NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))-pi/2*cos(pi/2*(x(:,1)-x(:,2))) ...  % Right hand side source
                             -NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))+pi/2*cos(pi/2*(x(:,1)-x(:,2)))];
    GD_Handle = @(x,varargin)[cos(pi/2*(x(:,1)+x(:,2))) -cos(pi/2*(x(:,1)+x(:,2)))];
    UEX_Handle1 = @(x,varargin)cos(pi/2*(x(:,1)+x(:,2)));
    UEX_Handle2 = @(x,varargin)-cos(pi/2*(x(:,1)+x(:,2)));
    UEX_GRAD_Handle1 = @(x,varargin)[-pi/2*sin(pi/2*(x(:,1)+x(:,2))) -pi/2*sin(pi/2*(x(:,1)+x(:,2)))];
    UEX_GRAD_Handle2 = @(x,varargin)[pi/2*sin(pi/2*(x(:,1)+x(:,2))) pi/2*sin(pi/2*(x(:,1)+x(:,2)))];
    UEX_Handle_P = @(x,varargin)sin(pi/2*((x(:,1)-x(:,2))));
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2_U = zeros(NREFS,1);
    Err_HCS_U = zeros(NREFS,1);
    Err_L2_P = zeros(NREFS,1);
    
    % Initialize mesh
    
    R_Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
    R_Mesh.ElemFlag = ones(size(R_Mesh.Elements,1),1);
    R_Mesh = add_Edges(R_Mesh);
    Loc = get_BdEdges(R_Mesh);
    R_Mesh.BdFlags = zeros(size(R_Mesh.Edges,1),1);
    R_Mesh.BdFlags(Loc) = -1;
    R_Mesh = refine_REG(R_Mesh);
    
    for i=1:NREFS
        
        R_Mesh = refine_REG(R_Mesh);        
        Loc = get_BdEdges(R_Mesh);
        Loc = unique([R_Mesh.Edges(Loc,1); R_Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(R_Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        Mesh = jiggle(R_Mesh,FixedPos);
        
        M_W(i) = get_MeshWidth(Mesh);

        % Solve the system

        A = assemMat_Stokes_CRP0(Mesh,@STIMA_Stokes_CRP0,NU,P7O6());
        L = assemLoad_Stokes_CRP0(Mesh,P7O6(),F_Handle);
        [U,FreeDofs] = assemDir_Stokes_CRP0(Mesh,-1,GD_Handle);
        L = L - A*U;
        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);

        % Extract each part from the solution vector
        
        nCoordinates = size(Mesh.Coordinates,1);
        nEdges = size(Mesh.Edges,1);
        nElements = size(Mesh.Elements,1);
        U_V1 = U(1:nEdges);
        U_V2 = U(nEdges+(1:nEdges));
        U_P = U(2*nEdges+(1:nElements));
        
        % Compute discretized error
        
        L2_U1 = L2Err_CR(Mesh,U_V1,P7O6(),UEX_Handle1);
        L2_U2 = L2Err_CR(Mesh,U_V2,P7O6(),UEX_Handle2);
        Err_L2_U(i) = sqrt(L2_U1^2+L2_U2^2);
        HCS_U1 = H1SErr_CR(Mesh,U_V1,P7O6(),UEX_GRAD_Handle1); 
        HCS_U2 = H1SErr_CR(Mesh,U_V2,P7O6(),UEX_GRAD_Handle2); 
        Err_HCS_U(i) = sqrt(HCS_U1^2+HCS_U2^2);
        Err_L2_P(i) = L2Err_P0(Mesh,U_P,P7O6(),UEX_Handle_P);
        save rate_CRP0_NC_Sqr.mat M_W Err_L2_U Err_HCS_U Err_L2_P
        
    end
    
    figure
    
    plot(M_W,Err_L2_U,'-^',M_W,Err_HCS_U,'-*',M_W,Err_L2_P,'-+');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p1 = polyfit(log(M_W),log(abs(Err_L2_U)),1);
    p2 = polyfit(log(M_W),log(abs(Err_L2_P)),1);
    add_Slope(gca,'SouthEast',p1(1));
    add_Slope(gca,'SouthWest',p2(1));
    ylabel('\bf\fontsize{14}Discretization error [log]')
    xlabel('\bf\fontsize{14}Mesh width [log]')
    title('\bfConvergence rate of discretization error for Crouzeix-Raviart elements')
    legend('\bf\fontsize{10}V, L2 Norm','\bf\fontsize{10}V, H1 Semi-norm',...
           '\bf\fontsize{10}P, L2 Norm');
    
    % Generate .eps files
    
    print('-depsc', 'rate_CRP0_NC_Sqr.eps');
    
    % Clear memory
    
%     clear all;