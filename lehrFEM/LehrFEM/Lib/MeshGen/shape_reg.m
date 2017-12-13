function varargout = shape_reg(Mesh)
% SHAPE_REG shape regularity
%
%   SR = SHAPE_REG(MESH) computes the value of the shape regularity
%   constant for a given mesh.
%
%   [SR,SR_DIST] = SHAPE_REG(MESH) also returns the values of the element
%   shape regularity constants of the mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of
%                the mesh.
%
%   Example:
%
%   sr = shape_reg(Mesh);

%   Copyright 2005-2005 Patrick Meury and Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Extract vetices
    
    vidx = Mesh.Elements;                   
    X1 = Mesh.Coordinates(vidx(:,1),1);     
    X2 = Mesh.Coordinates(vidx(:,2),1);
    X3 = Mesh.Coordinates(vidx(:,3),1);
    Y1 = Mesh.Coordinates(vidx(:,1),2);
    Y2 = Mesh.Coordinates(vidx(:,2),2);
    Y3 = Mesh.Coordinates(vidx(:,3),2);
    
    % Compute edge length
    
    E1 = sqrt((X2-X3).^2+(Y2-Y3).^2);      
    E2 = sqrt((X3-X1).^2+(Y3-Y1).^2);
    E3 = sqrt((X1-X2).^2+(Y1-Y2).^2);
    Radius = ((X1-X3).*(Y2-Y3)-(X2-X3).*(Y1-Y3))./(E1+E2+E3);   
    
    % Compute element shape regularity measure
    
    SR_Dist = max([E1 E2 E3],[],2)./Radius;
    
    % Compute shape regularity measure
    
    SR = max(SR_Dist);
    
    % Assign output arguement
    
    if(nargout == 1)  
      varargout{1} = SR;  
    elseif(nargout == 2)
      varargout{1} = SR;
      varargout{2} = SR_Dist;
    else
      error('Too many output arguments.');      
    end
    
return