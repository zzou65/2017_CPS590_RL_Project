% Run script for the convergence rate for strong regularization LFE2 finite
% element on the unit disk domain.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear Mesh;
    
    NREFS = 6;
    DHANDLE = @dist_circ;       % Signed distance function for a circle
    Center = [0 0];             % Center of circle
    R = 1;                      % Radius of circle

    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)pi^2*[cos(pi*x(:,2)) cos(pi*x(:,1))];
    GD_Handle = @(x,varargin)[cos(pi*x(:,2)) cos(pi*x(:,1))];
    UEX_Handle = @(x,varargin)[cos(pi*x(:,2)) cos(pi*x(:,1))];
    UEX_CURL_Handle = @(x,varargin)-pi*sin(pi*x(:,2))+pi*sin(pi*x(:,1));
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,1);
    Err_HCS = zeros(NREFS,1);
    
    % Initialize mesh
    
    Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    BdFlags = -1;
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = BdFlags;
    Mesh = refine_REG(Mesh,DHANDLE,Center,R);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh,DHANDLE,Center,R);
        
        M_W(i) = get_MeshWidth(Mesh);

        % Assemble Curl-curl matrix, MASS matrix and load vector

        nCoordinates = size(Mesh.Coordinates,1);
        [IC,JC,C] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
        [ID,JD,D] = assemMat_LFE2(Mesh,@STIMA_Div_LFE2,U_Handle,P7O6());
        [U,g,FreeDofs,IB,JB,B] = assemDir_StrRegLFE2(Mesh,BdFlags,GD_Handle);
        A = sparse([IC;ID;IB+2*nCoordinates;JB],[JC;JD;JB;IB+2*nCoordinates],[C;D;B;B]);
        l = assemLoad_LFE2(Mesh,P7O6(),F_Handle);
        L = [l;g];
        L = L - A*U;

        % Solve the system

        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        U = U(1:2*nCoordinates);

        % Compute discretization error
        
        Err_L2(i) = L2Err_LFE2(Mesh,U,P7O6(),UEX_Handle);
        Err_HCS(i) = HCurlSErr_LFE2(Mesh,U,P7O6(),UEX_CURL_Handle);

    end
    
    figure
    
    plot(M_W,Err_L2,'-^',M_W,Err_HCS,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_HCS)),1);
    add_Slope(gca,'West',p(1));
    ylabel('\bf\fontsize{14}Discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate of discretization error for strong regularization on the disk domain')
    legend('\bf\fontsize{10}StrReg:L2 Norm','\bf\fontsize{10}StrReg:H(Curl) Semi-Norm')
    
    % Generate .eps files
    
    print('-depsc', 'rate_StrReg_LFE2_Disk.eps');
    !gv rate_StrReg_LFE2_Disk.eps &
    
    % Clear memory
    
    clear all;
    
    