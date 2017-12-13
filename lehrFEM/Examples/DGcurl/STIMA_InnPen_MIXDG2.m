function Jloc = STIMA_InnPen_MIXDG2(Edge,Normal,BdFlag,LData,RData,varargin)
% STIMA_INNPEN_MIXDG2 Element stiffness matrix for interior penalty term.
%
%   JLOC = STIMA_INNPEN_MIXDG2(EDGE,NORMAL,BDFLAG,LDATA,RDATA,SHANDLE)
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
%   JLOC = STIMA_INNPEN_DGLFE(EDGE,NORMAL,BDFLAG,LDATA,RDATA,SHANDLE,SPARAM)
%   also handles the variable length argumet list SPARAM to the function
%   pointer SHANDLE.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Jloc = zeros(12,6);
  Jtmp = zeros(6,6);
  
  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);

  % Compute values of shape functions
     
  dS = norm(P1-P0);
  tmp = [1/3 1/6 0; 1/6 1/3 0; 0 0 0]/2;
  
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
  
  Jtmp = [ 1*tmp(LPerm,LPerm)  -1*tmp(LPerm,RPerm);
           1*tmp(RPerm,LPerm)  -1*tmp(RPerm,RPerm)]*dS;
            
  for j = 1:6 
      for i = 1:6
        Jloc(2*i-1:2*i,j) = Jtmp(i,j);
      end
  end
  
  for i = 1:2:12
        Jloc(i,:) = Normal(1)*Jloc(i,:);
        Jloc(i+1,:) = Normal(2)*Jloc(i+1,:);
  end
      
  if(LData.Element < RData.Element)
    Jloc = -Jloc;
  end
  
return
