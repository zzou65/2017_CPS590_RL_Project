function L = assemDir_dual(Mesh,BdFlags,L,FHandle,varargin)
% ASSEMDIR_dual Direchlet boundary conditions.
%
%   L = ASSEMDIR_DUAL(MESH,BDFLAG,L,QUADRULE,FHANDLE) incoporates the
%   Direchlet boundary conditions with the data given by FHANDLE into the
%   right hand side load vector L for the dual laplace problem. The boundary 
%   condition is only enforced
%   at edges whose BdFlag is equal to the integer BDFLAG. The 1D struct
%   QUADRULE is used to do the numerical integration along the edges.
%
%   L = ASSEMDIR_DUAL(MESH,BDFLAG,L,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the boundary data function
%   FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first column
%                 specifies the left hand side element where the second column
%                 specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of elements. 
%
%   Example:
%
%   L = assemDir_Dual(Mesh,BdFlags,L,QuadRule,FHandle);
%
%   Copyright 2005-2006 Patrick Meury & Kah-Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  Lloc = zeros(1,1);
  
    Loc = get_BdEdges(Mesh);
    for i = Loc'
        vid=Mesh.Edges(i,:);
        
        if(Mesh.Edge2Elem(i,1)>0)
           p=1;
        else p=-1;
        end
        a1=Mesh.Coordinates(vid(1),:);
        a2=Mesh.Coordinates(vid(2),:);
        l=norm(a2-a1);
        
        QuadRule=gauleg(0,1,10);
        x=ones(size(QuadRule.x,1),1)*a1+QuadRule.x*(a2-a1);
              
        % Evaluate Direchlet boundary data
      
        FVal = FHandle(x,-1,varargin{:});
        
       % Numerical integration along an edge
    
        Lloc = sum(QuadRule.w.*FVal*p);
        
      % Add contributions of Neumann data to load vector
      
      L(i) = L(i)+Lloc;
      
    end    
  
return