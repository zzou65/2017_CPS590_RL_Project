function Jloc = STIMA_BndPen_DGCR(Edge,Normal,BdFlag,Data,SHandle,varargin)
% STIMA_BNDPEN_DGCR Element stiffness matrix for boundary penalty term.
%
%   JLOC = STIMA_BNDPEN_DGCR(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,SHANDLE)
%   computes the entries of the element stiffness matrix for the boundary
%   penalty term.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   The integer BDFLAG denotes the boundary flag of the current edge. Note
%   that for interior edges only values larger than are allowed.
%
%   The struct DATA contains the left or right hand side element data:
%    ELEMENT  Integer specifying the neighbouring element.
%    ELEMFLAG Integer specifying the element flag of the neighbouring
%             element or zero.
%    VERTICES 3-by-2 matrix specifying the vertices of the neighbouring
%             element.
%    EDGELOC  Integer specifying the local edge number on the neighbouring
%             element.
%    MATCH    Integer specifying the relative orientation of the edge with
%             respect to the orientation of the neighbouring element.
%
%   SHANDLE is a function pointer to the edge weight function.
%
%   JLOC = STIMA_BNDPEN_DGCR(EDGE,NORMAL,BDFLAG,LDATA,RDATA,SHANDLE,SPARAM)
%   also handles the variable length argumet list SPARAM to the function
%   pointer SHANDLE.
%
%   Example:
%
%   Sigma = @(P1,P0,varargin)norm(P1-P0);
%   Jloc = STIMA_BndPen_DGCR(Edge,Normal,BdFlags,Data,Sigma);
%
%   See also shap_DGCR.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Jloc = zeros(3,3);

  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  sigma = SHandle(P0,P1,varargin{:});
  
  % Compute values of shape functions
    
  dS = norm(P1-P0);
  tmp = [1/3 -1/3 0; -1/3 1/3 0; 0 0 1];
  
  switch(Data.EdgeLoc)
    case 1
      if(Data.Match == 1)
        Perm = [3 1 2];  
      else
        Perm = [3 2 1];  
      end
    case 2
      if(Data.Match == 1)
        Perm = [2 3 1];  
      else
        Perm = [1 3 2];  
      end
    case 3 
      if(Data.Match == 1)
        Perm = [1 2 3];  
      else
        Perm = [2 1 3];  
      end
  end
  
  Jloc = sigma*tmp(Perm,Perm)*dS;
      
return
