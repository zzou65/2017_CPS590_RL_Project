function Mesh = orient_Elems(Mesh)
% ORIENT_ELEMS check and correct the orientation of each element 
%
%   MESH = ORIENT_ELEMS(MESH) check the orientation of each element and 
%   change the clock-wise one to counter-clock-wise
%
%   Example:
%
%   Mesh = orient_Elems(Mesh);

%   Copyright 2005-2007 Patrick Meury & Mengyu Wang & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Change orientation of all elements to counter-clockwise
  
  for i = 1:size(Mesh.Elements,1);
    v = Mesh.Coordinates(Mesh.Elements(i,2),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    w = Mesh.Coordinates(Mesh.Elements(i,end),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    if(v(1)*w(2)-v(2)*w(1) < 0)
      Mesh.Elements(i,:) = Mesh.Elements(i,end:-1:1);
    end
  end  