function arg_out = rotate(arg_in,phi)
% ROTATE Rotates coordinates.
%
%   COORDINATES = ROTATE(COORDINATES,PHI) rotates all vertices in the array
%   COORDINATES by the angle PHI in a counter-clockwise orientation.
%
%   MESH = ROTATE(MESH,PHI) rotates all vertices in the struct MESH by the
%   angle PHI in a counter-clockwise orientation.
%
%   Example:
%
%   Mesh = rotate(Mesh,3/2*pi);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Set up rotation matrix
  
  A = [cos(phi) sin(phi); -sin(phi) cos(phi)];
 
  % Rotate coordinates and assign output arguments
  
  arg_out = deal(arg_in);
  if(isstruct(arg_in))
    arg_out.Coordinates = arg_in.Coordinates*A;    
  else
    arg_out = arg_in*A;
  end
  
return