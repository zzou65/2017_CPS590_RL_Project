function arg_out = stretch(arg_in,x_dir,y_dir)
% STRETCH Stretches coordinates.
%
%   COORDINATES = STRETCH(COORDINATES,X_DIR,Y_DIR) stretches all vertices in
%   the array COORDINATES by X_DIR in x-direction and by Y_DIR in y-direction.
%
%   MESH = STRETCH(MESH,X_DIR,Y_DIR) stretches all vertices in the struct MESH
%   by X_DIR in x-direction and Y_DIR in y-Direction.
%
%   Example:
%
%   Mesh = stretch(Mesh,2,3);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Stretch all coordinates
  
  arg_out = deal(arg_in);
  if(isstruct(arg_in))
    arg_out.Coordinates(:,1) = x_dir*arg_in.Coordinates(:,1);
    arg_out.Coordinates(:,2) = y_dir*arg_in.Coordinates(:,2);
  else
    arg_out(:,1) = x_dir*arg_in(:,1);
    arg_out(:,2) = y_dir*arg_in(:,2);
  end
  
return