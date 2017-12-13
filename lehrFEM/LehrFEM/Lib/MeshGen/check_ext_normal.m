function orient = check_ext_normal(Normal,EdgeVert,ElemVert,EdgeLoc)
%CHECK_EXT_NORMAL check orientation of normal vector
%
%   ORIENT = CHECK_EXT_NORMAL(NORMAL,EDGEVERT,ELEMVERT,EDGELOC) checks if
%   NORMAL is an exterior normal vector of the edge EDGEVERT of the element
%   given by ELEMVERT.  ORIENT is 1 if this is the case and -1 otherwise,
%   so ORIENT*NORMAL is an exterior normal vector.
%
%   EDGEVERT is a 2-by-2 matrix containing the coordinates of the endpoints
%   of the edge in question in its rows.
%
%   ELEMVERT is a 3-by-2 matrix containing the coordinates of the vertices
%   of the triangular element in its rows.
%
%   EDGELOC is the row index of the vertex of the triangle not on the
%   current edge.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%   if(~(isequal(ElemVert(setdiff(1:3,EdgeLoc),:),EdgeVert)||isequal(ElemVert(setdiff(1:3,EdgeLoc),:),EdgeVert([2 1],:))))
%     error('edge and element do not match');
%   end

  intPnt = ElemVert(EdgeLoc,:);         % point inside element
  outVec = EdgeVert(1,:) - intPnt;      % vector crossing edge from inside to outside
  orient = sign(outVec*Normal');        % does the normal vector cross the edge in the same direction?

return