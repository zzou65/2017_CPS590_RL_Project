function mean = mean_LFE(U,Mesh)
% MEAN_LFE computes the mean for linear finite elements.
%
%   MEAN = MEAN_LFE(U,MESH) computes the mean values for the linear 
%   finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   mean = mean_LFE(U,MESH);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nElements = size(Mesh.Elements,1);
    
  % Compute mean value

  mean = 0;
  omega = 0;
  
  for i = 1:nElements
       
    % Extract vertex numbers
    
    i1 = Mesh.Elements(i,1);
    i2 = Mesh.Elements(i,2);
    i3 = Mesh.Elements(i,3);
      
    % Compute element mapping  
        
    BK = [Mesh.Coordinates(i2,:)-Mesh.Coordinates(i1,:); ...
          Mesh.Coordinates(i3,:)-Mesh.Coordinates(i1,:)];
    det_BK = abs(det(BK));
      
    % Compute mean at current element 
      
    mean = mean + sum(U(Mesh.Elements(i,:)))*det_BK/6;
    omega = omega + det_BK/2;
    
  end
 
  mean = mean/omega;
  
return