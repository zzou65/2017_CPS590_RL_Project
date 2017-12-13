function varargout = limit_ad(Mesh,u_old,alpha,varargin)
% LIMIT_AD Adapative, piecewise constant limiter.
%
%   U = LIMIT_P0(MESH,U,ALPHA) projects the piecewie linear function U onto
%   the space of piecewise constant functions.
%
%   ALPHA denotes the exponent of the volume weight inside the indicator
%   function and must be between 2 and 6.
%
%   Example:
%
%   U = limit_ad(Mesh,U,3);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Preallocate memory
  
  u_new = zeros(3*nElements,1);
  
  % Extract boundary flags
  
  BdFlags = Mesh.BdFlags;
  
  % Compute limited value
  
  beta = (alpha-2)/3;
  u_l = zeros(1,3);
  u_r = zeros(1,3);
  gK = zeros(nElements,1);
  for i = 1:nElements
           
    % Extract current element data 
      
    vidx = Mesh.Elements(i,:);
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
          Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK)); 
      
    for j = 1:3
        
      % Extract current edge data
      
      id_s = Mesh.Elements(i,rem(j,3)+1);
      id_e = Mesh.Elements(i,rem(j+1,3)+1);
      eidx = Mesh.Vert2Edge(id_s,id_e);
                  
      if(BdFlags(eidx) >= 0)
      
        % Extract left and right hand side neighbouring element data  
         
        ElemLeft = Mesh.Edge2Elem(eidx,1);
        LeftLoc = Mesh.EdgeLoc(eidx,1);
        
        ElemRight = Mesh.Edge2Elem(eidx,2);
        RightLoc = Mesh.EdgeLoc(eidx,2); 
                 
        % Compute indicator value on the current element
      
        idx_l = 3*(ElemLeft-1)+[1 2 3];
        switch(LeftLoc)
          case 1
            if(Mesh.EdgeOrient(eidx,1) == 1)
              u_l(3) = u_old(idx_l(1));
              u_l(1) = u_old(idx_l(2));
              u_l(2) = u_old(idx_l(3));
            else
              u_l(3) = u_old(idx_l(1));
              u_l(2) = u_old(idx_l(2));
              u_l(1) = u_old(idx_l(3));
            end
          case 2
            if(Mesh.EdgeOrient(eidx,1) == 1)
              u_l(2) = u_old(idx_l(1));
              u_l(3) = u_old(idx_l(2));
              u_l(1) = u_old(idx_l(3));
            else
              u_l(1) = u_old(idx_l(1));
              u_l(3) = u_old(idx_l(2));
              u_l(2) = u_old(idx_l(3));  
            end
          case 3
            if(Mesh.EdgeOrient(eidx,1) == 1)
              u_l(1) = u_old(idx_l(1));
              u_l(2) = u_old(idx_l(2));
              u_l(3) = u_old(idx_l(3));
            else
              u_l(2) = u_old(idx_l(1));
              u_l(1) = u_old(idx_l(2));
              u_l(3) = u_old(idx_l(3));
            end
        end
      
        idx_r = 3*(ElemRight-1)+[1 2 3];
        switch(RightLoc)
          case 1
            if(Mesh.EdgeOrient(eidx,2) == 1)
              u_r(3) = u_old(idx_r(1));
              u_r(1) = u_old(idx_r(2));
              u_r(2) = u_old(idx_r(3));
            else
              u_r(3) = u_old(idx_r(1));
              u_r(2) = u_old(idx_r(2));
              u_r(1) = u_old(idx_r(3));  
            end
          case 2
            if(Mesh.EdgeOrient(eidx,2) == 1)
              u_r(2) = u_old(idx_r(1));
              u_r(3) = u_old(idx_r(2));
              u_r(1) = u_old(idx_r(3));
            else
              u_r(1) = u_old(idx_r(1));
              u_r(3) = u_old(idx_r(2));
              u_r(2) = u_old(idx_r(3));  
            end
          case 3
            if(Mesh.EdgeOrient(eidx,2) == 1)
              u_r(1) = u_old(idx_r(1));
              u_r(2) = u_old(idx_r(2));
              u_r(3) = u_old(idx_r(3));
            else
              u_r(2) = u_old(idx_r(1));
              u_r(1) = u_old(idx_r(2));
              u_r(3) = u_old(idx_r(3));  
            end
        end

        v = u_l-u_r;
                
        gK(i) = gK(i) + 1/3*(v(1)-v(2))^2+v(3)^2;
      
      end
            
    end
      
    gK(i) = (2/det_BK)^beta*gK(i);
    
    % Check indicator criterion
    
    idx = 3*(i-1)+[1 2 3];
    if(gK(i) > 1)       
      u_new(idx) = sum(u_old(idx))/3; 
    else
      u_new(idx) = u_old(idx);  
    end
       
  end
    
  % Assign output arguments
  
  varargout{1} = u_new;
  if(nargout > 1)
    varargout{2} = gK; 
  end
  
return