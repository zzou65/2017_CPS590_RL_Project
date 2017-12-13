function arg_out = shift(arg_in,x0)
% SHIFT Shifts coordinates.
%
%   COORDINATES = SHIFT(COORDINATES,X0) shifts all vertices in the array
%   COORDINATES by the vector X0.
%
%   MESH = SHIFT(MESH,X0) shifts all vertices in the struct MESH by
%   the vector X0.
%
%   Example:
%
%   Mesh = shift(Mesh,[1 2]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Shift all coordinates
  
  arg_out = deal(arg_in);
  if(isstruct(arg_in))
    arg_out.Coordinates(:,1) = arg_in.Coordinates(:,1)+x0(1);
    arg_out.Coordinates(:,2) = arg_in.Coordinates(:,2)+x0(2);
  else
    arg_out(:,1) = arg_in(:,1)+x0(1);
    arg_out(:,2) = arg_in(:,2)+x0(2);
  end
  
return