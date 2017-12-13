function Jloc = STIMA_InnPen_DGCR(Edge,Normal,BdFlag,LData,RData,SHandle,varargin)
% STIMA_INNPEN_DGCR Element stiffness matrix for interior penalty term.
%
%   JLOC = STIMA_INNPEN_DGCR(EDGE,NORMAL,BDFLAG,LDATA,RDATA,SHANDLE)
%   computes the entries of the element stiffness matrix for the interior
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
%   The structs LDATA and RDATA conatin the left and right hand side
%   element data:
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
%   JLOC = STIMA_INNPEN_DGCR(EDGE,NORMAL,BDFLAG,LDATA,RDATA,SHANDLE,SPARAM)
%   also handles the variable length argumet list SPARAM to the function
%   pointer SHANDLE.
%
%   Example:
%
%   Sigma = @(P1,P0,varargin)norm(P1-P0);
%   Jloc = STIMA_InnPen_DGCR(Edge,Normal,BdFlags,LDatat,RData,Sigma);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Jloc = zeros(6,6);
  
  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  sigma = SHandle(P0,P1,varargin{:});

  % Compute values of shape functions
 
  dS = norm(P1-P0);
  tmp = [1/3 -1/3 0; -1/3 1/3 0; 0 0 1];
  
  switch(LData.EdgeLoc)
    case 1
      if(LData.Match == 1)
        LPerm = [3 1 2];  
      else
        LPerm = [3 2 1];  
      end
    case 2
      if(LData.Match == 1)
        LPerm = [2 3 1];  
      else
        LPerm = [1 3 2];  
      end
    case 3
      if(LData.Match == 1)
        LPerm = [1 2 3];  
      else
        LPerm = [2 1 3];  
      end
  end
  
  switch(RData.EdgeLoc)
    case 1
      if(RData.Match == 1)
        RPerm = [3 1 2];  
      else
        RPerm = [3 2 1];  
      end
    case 2
      if(RData.Match == 1)
        RPerm = [2 3 1];
      else
        RPerm = [1 3 2];  
      end
    case 3
      if(RData.Match == 1)
        RPerm = [1 2 3];  
      else
        RPerm = [2 1 3];  
      end
  end
  
  Jloc = sigma*[ tmp(LPerm,LPerm) -tmp(LPerm,RPerm); ...
                -tmp(RPerm,LPerm)  tmp(RPerm,RPerm)]*dS;
    
return
