function varargout = AssemDir_StrRegLFE2(Mesh,BdFlags,GD_Handle,varargin)
% ASSEMDIR_LFE2 Dirichlet boundary conditions for strong regularization 
%
%   [U,G,FREEDOFS,B] = ASSEMDIR_STRREGLFE2(MESH,BDFLAG,GD_HANDLE) 
%   incoporates the Dirichlet boundary conditions with the data given by 
%   GD_HANDLE into the finite element solution U. The boundary condition is
%   only enforced at the corner vertices.
%   G is the Dirichelet boundary data.
%   FreeDofs is all the non-corner node as degree of freedom.
%   B is the matrix for constraints of non-corner boundary Dofs.
%
%   [U,G,FREEDOFS,B] = ASSEMDIR_STRREGLFE2(MESH,BDFLAG,GD_HANDLE,EPARAM)
%   handles the variable length argument list EPARAM to the function handle
%   GD_HANDLE during the assembly process.
%
%   [U,G,FREEDOFS,IB,JB,KB] = ASSEMDIR_STRREGLFE2(MESH,BDFLAG,GD_HANDLE)
%   assembles the constraint matrix B in an array representation.
%
%   [U,G,FREEDOFS,B,NCDOFS] = ASSEMDIR_STRREGLFE2(MESH,BDFLAG,GD_HANDLE,EPARAM)
%   returns non-cornered D.O.Fs.
%
%   [U,G,FREEDOFS,IB,JB,KB,NCDOFS] = ASSEMDIR_STRREGLFE2(MESH,BDFLAG,...
%   GD_HANDLE,EPARAM) returns B in array representation and also returns
%   non-cornered D.O.Fs.   
% 
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%
%   Example:
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
%   Mesh = add_Edges(Mesh);
%   Mesh = add_Edge2Elem(Mesh);
%   Loc = get_BdEdges(Mesh);
%   BdFlags = [-1 -2 -3 -4 -5 -6];
%   Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
%   Mesh.BdFlags(Loc)= BdFlags;
%   Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
%   GD_Handle = @(x,varargin)[cos(pi*x(:,2)) cos(pi*x(:,1))];
%   [U,g,FreeDofs,IB,JB,B] = assemDir_StrRegLFE2(Mesh,BdFlags,GD_Handle);
%  
%   See also ASSEMMAT_LFE2.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    nCoordinates = size(Mesh.Coordinates,1);
    
    % Preallocate memory
    
    TanInfo = zeros(nCoordinates,2);
    FlagInfo = zeros(nCoordinates,2);
    u = zeros(2*nCoordinates,1);
    g = zeros(nCoordinates,1);
    B = zeros(2*nCoordinates,1);
    CDofs = [];
    
    % Extract doundary information from the mesh
    
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    BdNodes = unique(Mesh.Edges(Loc,:))';
    
    for j = BdFlags
        
        DEdges = Loc(Mesh.BdFlags(Loc) == j)';
        
        for i = DEdges

            % Extract edge nodes
            vidx1 = Mesh.Edges(i,1); vidx2 = Mesh.Edges(i,2);
            % Match orientation
            if (Mesh.Edge2Elem(i,1) ~= 0 & Mesh.Edge2Elem(i,2) == 0)
                sign = 1;
            else 
                sign = -1;
            end            
            % Record information            
            tan = (Mesh.Coordinates(vidx2,:) - Mesh.Coordinates(vidx1,:))*sign*.5;
            TanInfo(vidx1,:) = TanInfo(vidx1,:) + tan;
            TanInfo(vidx2,:) = TanInfo(vidx2,:) + tan;            
            if ( FlagInfo(vidx1,1) == 0 )
                FlagInfo(vidx1,1) = j;                
            else                
                FlagInfo(vidx1,2) = j;                
            end            
            if ( FlagInfo(vidx2,1) == 0 )                
                FlagInfo(vidx2,1) = j;                
            else                
                FlagInfo(vidx2,2) = j;                
            end
            
        end
        
    end
    
    % Assemble boundary data
    
    I = [1:nCoordinates;1:nCoordinates];
    I = I(:);
    J = [1:nCoordinates;nCoordinates+1:2*nCoordinates];
    J = J(:);
   
    for i = BdNodes
        
        GVal = GD_Handle(Mesh.Coordinates(i,:));
        
        if (FlagInfo(i,1)==FlagInfo(i,2))
            B([2*i-1 2*i]) = TanInfo(i,:);
            g(i) = dot(GVal,TanInfo(i,:));
        else
            CDofs = [CDofs;i];
            u([i i+nCoordinates]) = GVal; 
        end
        
    end
    
    NCNodes = setdiff(1:nCoordinates,CDofs)'; % Non-cornered nodes
    NCDofs = [NCNodes;NCNodes+nCoordinates];  % Non-cornered D.O.Fs
    FDofs = setdiff(BdNodes,CDofs);           % Boundary but non-corner nodes
    B = sparse(I,J,B);
    B = B(FDofs,:);
    g = g(FDofs);
    [I,J,B] = find(B);
    nFDofs = length(FDofs);
    FreeDofs = [NCNodes;NCNodes+nCoordinates;(2*nCoordinates+1:2*nCoordinates+nFDofs)'];
    U = zeros(2*nCoordinates+nFDofs,1);
    U(1:2*nCoordinates) = u;
    
    % Assign output arguments

    if(nargout == 6)
        varargout{1} = U;
        varargout{2} = g;
        varargout{3} = FreeDofs;
        varargout{4} = I;
        varargout{5} = J;
        varargout{6} = B;
    elseif(nargout == 4)
        varargout{1} = U;
        varargout{2} = g;
        varargout{3} = FreeDofs;
        varargout{4} = sparse(I,J,B);
    elseif(nargout == 7)
        varargout{1} = U;
        varargout{2} = g;
        varargout{3} = FreeDofs;
        varargout{4} = I;
        varargout{5} = J;
        varargout{6} = B;
        varargout{7} = NCDofs;
    elseif(nargout == 5)
        varargout{1} = U;
        varargout{2} = g;
        varargout{3} = FreeDofs;
        varargout{4} = sparse(I,J,B);
        varargout{5} = NCDofs;
        
    end
  
return