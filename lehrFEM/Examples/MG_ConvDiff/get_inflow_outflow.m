function [inflow,outflow,neutral,v] = get_inflow_outflow(Mesh,v_handle)
%GET_INFLOW_OUTFLOW get inflow and outlfow edges 
%   
%   [INFLOW,OUTFLOW,NEUTRAL] = GET_INFLOW_OUTFLOW(MESH,V_HANDLE)
%   determines the indices of the inflow, outflow and neutral edges of the
%   velocity field V_HANDLE on the meeh MESH, stored, respectively, in the
%   vectors INFLOW, OUTFLOW and NEUTRAL.
%
%   [INFLOW,OUTFLOW,NEUTRAL,VEL] = GET_INFLOW_OUTFLOW(...) also returns the
%   velocities VEL at the midpoints of the edges.
%
%   Example : 
%
%   Mesh = add_Edges(load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'));
%   [inflow,outflow,neutral] = get_inflow_outflow(Mesh,@(x,varargin) ones(size(x)));

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % find boundary edges

  Loc = get_BdEdges(Mesh);
  
  % determine velocity on each edge
  
  midpt = 0.5*(Mesh.Coordinates(Mesh.Edges(Loc,1),:)+Mesh.Coordinates(Mesh.Edges(Loc,2),:));
  v = v_handle(midpt);
  
  % generate vector that points outside
  
  if(~isfield(Mesh,'EdgeLoc'))
    Mesh = add_Edge2Elem(Mesh);
  end
  EdgeLoc = max(Mesh.EdgeLoc(Loc,:),[],2);
  Edge2Elem = max(Mesh.Edge2Elem(Loc,:),[],2);
  inpt = Mesh.Coordinates(diag(Mesh.Elements(Edge2Elem,EdgeLoc)),:);
  outvec = midpt-inpt;
  
  % calculate some normal vectors
  
  edge = Mesh.Coordinates(Mesh.Edges(Loc,2),:)-Mesh.Coordinates(Mesh.Edges(Loc,1),:);
  n = [-edge(:,2),edge(:,1)];
  
  % make normals point outward
  
  s = sign(sum(outvec.*n,2));
  n = s(:,[1 1]).*n;
  
  % calculate scalar products with velocities
  
  scapro = sum(v.*n,2);
  
  % determine inflow, outflow and neutral boundaries
  
  inflow = Loc(scapro<0);
  outflow = Loc(scapro>0);
  neutral = Loc(scapro==0);

return