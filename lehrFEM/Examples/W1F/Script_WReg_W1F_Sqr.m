% Run script for the convergence rate for W1F finite element solver on the
% unit square domain.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    NREFS = 3;
    
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_CURL_Handle = @(x,varargin)pi*cos(pi*x(:,2))-pi*cos(pi*x(:,1));
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,1);
    Err_HCS = zeros(NREFS,1);
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1; -1 1];
    Mesh.Elements=[1 2 4;2 3 4];
    
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=-1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh);
        M_W(i) = get_MeshWidth(Mesh);
        
        % Assemble Curl-curl matrix, MASS matrix and load vector

        [IC,JC,C] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
        B = assemMat_WRegW1F(Mesh,@STIMA_WReg_W1F);
        D = assemMat_LFE(Mesh,@MASS_Lump_LFE);
        
        nCoordinates = size(Mesh.Coordinates,1);
        Loc = get_BdEdges(Mesh);
        DEdges = Loc(Mesh.BdFlags(Loc) == -1);
        DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
        FreeDofs = setdiff(1:nCoordinates,DNodes);
        B = B(:,FreeDofs);
        D = D(FreeDofs,FreeDofs);
    
        T = B*inv(D)*transpose(B);
        [IT,JT,T] = find(T);
        A = sparse([IC;IT],[JC;JT],[C;T]);
        L = assemLoad_W1F(Mesh,P7O6(),F_Handle);


        % Incorporate Dirichlet boundary data

        [U,FreeDofs] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
        L = L - A*U;
        
        
        % Solve the system

        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        Err_L2(i) = L2Err_W1F(Mesh,U,P7O6(),UEX_Handle);
        Err_HCS(i) = HCurlSErr_W1F(Mesh,U,P7O6(),UEX_CURL_Handle);
        
    end
    
    figure
    
    plot(M_W,Err_L2,'-^',M_W,Err_HCS,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2)),1);
    add_Slope(gca,'SouthEast',p(1));
    ylabel('\bf\fontsize{14}Discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for discretization error on the square domain')
    legend('\bf\fontsize{10}WRegW1F:L2 Norm','\bf\fontsize{10}WRegW1F:H(Curl) Semi-Norm')
    
    % Generate .eps files
    
    print('-depsc', 'rate_WRegW1F_Sqr.eps');
    !gv rate_WRegW1F_Sqr.eps &
    
    % Clear memory
    
%     clear all;
    