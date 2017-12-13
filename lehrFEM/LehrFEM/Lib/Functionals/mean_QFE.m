function mean = mean_QFE(U,Mesh)
% MEAN_QFE computes the mean for quadratic finite element.
%
%   MEAN = MEAN_QFE(U,MESH) computes the mean values for the quadratic 
%   finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   mean = mean_QFE(U,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
    
  % Compute mean value

  mean = 0;
  omega = 0;
  
  for i = 1:nElements
       
    % Extract vertex numbers
    
    i1 = Mesh.Elements(i,1);
    i2 = Mesh.Elements(i,2);
    i3 = Mesh.Elements(i,3);
    
    i4 = Mesh.Vert2Edge(i1,i2) + nCoordinates;
    i5 = Mesh.Vert2Edge(i2,i3) + nCoordinates;
    i6 = Mesh.Vert2Edge(i3,i1) + nCoordinates;
      
    % Compute element mapping  
        
    BK = [Mesh.Coordinates(i2,:)-Mesh.Coordinates(i1,:); ...
          Mesh.Coordinates(i3,:)-Mesh.Coordinates(i1,:)];
    det_BK = abs(det(BK));
      
    % Compute mean at current element 
      
    mean = mean + (U(i4)+U(i5)+U(i6))*det_BK/6;
    omega = omega + det_BK/2;
    
  end
  
  mean = mean/omega;
    
return
