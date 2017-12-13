% Run script for strong regularization LFE2 finite element solver on the 
% unit square domian.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 3;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)[5*ones(size(x,1),1) ones(size(x,1),1)];
    GD_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];

    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    BdFlags = [-1 -2 -3 -4];
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=BdFlags;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh);
    Mesh = refine_REG(Mesh);
    
    for i=1:NREFS

        Mesh = refine_REG(Mesh);

        % Solve the system

        nCoordinates = size(Mesh.Coordinates,1);
        [IC,JC,C] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
        [ID,JD,D] = assemMat_LFE2(Mesh,@STIMA_Div_LFE2,U_Handle,P7O6());
        [U,g,FreeDofs,IB,JB,B] = assemDir_StrRegLFE2(Mesh,BdFlags,GD_Handle);
        A = sparse([IC;ID;IB+2*nCoordinates;JB],[JC;JD;JB;IB+2*nCoordinates],[C;D;B;B]);
        l = assemLoad_LFE2(Mesh,P7O6(),F_Handle);
        L = [l;g];
        L = L - A*U;
        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        u = U(1:2*nCoordinates);

        % Plot the solution

        ux = u(1:nCoordinates);
        uy = u(nCoordinates+1:2*nCoordinates);
        norm_u = (ux.^2+uy.^2).^.5;

        fig = plot_LFE(norm_u,Mesh);
        set(gcf,'renderer','zbuffer')
        colorbar
        title(['{\fontsize{10}\bfNodal Elements, level' int2str(i) '}'])
        FileName = ['Distri_LFE2_Sqr_' int2str(i) '.eps'];
        print('-depsc', FileName);
        close(fig);
        system(['gv ' FileName ' &']);
        
    end

    
    % Clear memory
    
%     clear all
    
    