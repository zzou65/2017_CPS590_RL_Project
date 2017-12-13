% Run script for norm of the solution of LFE2 finite element on the unit
% disk domain.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 5;
    DHANDLE = @dist_circ;       % Signed distance function for a circle
    Center = [0 0];             % Center of circle
    R = 1;                      % Radius of circle
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)[5*ones(size(x,1),1) ones(size(x,1),1)];
    GD_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];
    Dummy_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];
    Dummy_C_Handle = @(x,varargin)zeros(size(x,1),1);
    
    % Preallocate memory

    Norm_L2_LFE2 = zeros(NREFS,1);
    Norm_HCS_LFE2 = zeros(NREFS,1);
    Norm_L2_W1F = zeros(NREFS,1);
    Norm_HCS_W1F = zeros(NREFS,1);
    
    % LFE2 results
    
    % Initialize mesh
    
    Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    BdFlags = -1;
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=BdFlags;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh,DHANDLE,Center,R);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh,DHANDLE,Center,R);

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
        u = U(1:2*nCoordinates);
        
        % Compute discretization error
        
        Norm_L2_LFE2(i) = L2Err_LFE2(Mesh,U,P7O6(),Dummy_Handle);
        Norm_HCS_LFE2(i) = HCurlSErr_LFE2(Mesh,U,P7O6(),Dummy_C_Handle);
        
    end
    
    % W1F results
       
    clear Mesh;
    
    % Initialize mesh
    
    Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=-1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh,DHANDLE,Center,R);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh,DHANDLE,Center,R);
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
        Norm_L2_W1F(i) = L2Err_W1F(Mesh,U,P7O6(),Dummy_Handle);
        Norm_HCS_W1F(i) = HCurlSErr_W1F(Mesh,U,P7O6(),Dummy_C_Handle);
        
    end
    
    % Plot the results
    
    plot(1:NREFS,Norm_L2_LFE2,'-^',1:NREFS,Norm_HCS_LFE2,'-*',...
         1:NREFS,Norm_L2_W1F,'->',1:NREFS,Norm_HCS_W1F,'-o');
    xlabel('\bf\fontsize{14}Level')
    title('\bfL^{2} norm of each solution on the unit disk domain')
    ylabel('\bf\fontsize{14}L^{2}-Norm')
    legend('\bf\fontsize{10}L^{2} nodal elements','\bf\fontsize{10}L^{2} curl nodal elements',...
           '\bf\fontsize{10}L^{2} edge elements','\bf\fontsize{10}L^{2} curl edge elements')
    axis([0,6,0,10])
    
    % Output .eps file
    
    print('-depsc','L2Norm_Disk.eps');
    !gv L2Norm_Disk.eps;
    
    % Clear memory
    
%     clear all