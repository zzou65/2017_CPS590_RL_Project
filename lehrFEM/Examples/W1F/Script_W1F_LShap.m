% Run script for the convergence rate for W1F finite element solver on the
% L-Shaped domain.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    NREFS = 2;
    
    MU_HANDLE=@(x,varargin)1;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @f_LShap_W1F;
    GD_Handle = @g_D_LShap_W1F;
    UEX_Handle = @uex_LShap_W1F;
    UEX_CURL_Handle = @uex_Curl_LShap;
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,1);
    Err_HCS = zeros(NREFS,1);
    
    
    Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
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
        [IM,JM,M] = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P3O3());
        A = sparse([IC;IM],[JC;JM],[C;M]);
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
    add_Slope(gca,'NorthWest',p(1));
    p = polyfit(log(M_W),log(abs(Err_HCS)),1);
    add_Slope(gca,'SouthEast',p(1));
    ylabel('\bf\fontsize{14}Discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for discretization error on the L-Shaped domain')
    legend('\bf\fontsize{10}W1F:L2 Norm','\bf\fontsize{10}W1F:H(Curl) Semi-Norm')
    
    % Generate .eps files
    
%     print('-depsc', 'rate_W1F_LShap.eps');
%     !gv rate_W1F_LShap.eps &
%     
    % Clear memory
    
    clear all;
    
    