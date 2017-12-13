function L = assemLoad_LFE_S1(Mesh,QuadRule,FHandle,varargin)
%ASSEMLOAD_LFE_S1 Assemble load vector for integral over unit circle
%   
%   L = ASSEMLOAD_LFE_S1(MESH,QUADRULE) assembles the global load vector
%   for the right hand side functional given by the integral of the trace
%   on the unit circle.
%
%    The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%   Note that QUADRULE must be a one dimensional quadrature rule on the
%   interval [0,1].
%
%   L = ASSEMLOAD_LFE_S1(MESH,QUADRULE,FHANDLE[,FPARAM]) uses FHANDLE as a
%   weight for the measure in the funcitonal.  FHANDLE must take the
%   arguments Phi,ElemFlag,FPARAM, where Phi is the standard
%   parametrization of the unit circle and ElemFlag specifies additional
%   element information.
%
%   Example:
%
%   QuadRule = gauleg(0,1,1);
%   L = assemLoad_LFE_S1(Mesh,QuadRule);

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Assign default values to omitted arguments
   
  if(nargin < 3 || isempty(FHandle))
    FHandle = @(x,varargin) ones(size(x));
  end

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    vcoord = Mesh.Coordinates(vidx,:);
    
    % Find points on intersection of unit circle and boundary of element
    
    endpts = zeros(0,2);
    
    for j = 1:3
      
      % Look for t such that ||v+tw|| = 1, 0<=t<=1.
      
      v = vcoord(mod(j-1,3)+1,:);
      w = vcoord(mod(j,3)+1,:)-v;
      
      a = w*w';
      b = 2*v*w';
      c = v*v'-1;
      
      if(a == 0)
        t = -c/b;
        if(t>=0 && t<=1)
          endpts = [endpts;v+t*w];
        end 
      else
        discr = b^2-4*a*c;

        if(discr == 0)
          t = -b/(2*a);
          if(t>=0 && t<=1)
            endpts = [endpts;v+t*w];
          end
        elseif(discr > 0)
          t = (-b-sqrt(discr))/(2*a);
          if(t>=0 && t<=1)
            endpts = [endpts;v+t*w];
          end
          t = (-b+sqrt(discr))/(2*a);
          if(t>=0 && t<=1)
            endpts = [endpts;v+t*w];
          end
        end
      end
      
    end % for loop over edges of triangle (index j)
    
    % Continue if there are no such points

    if(isempty(endpts))
      continue
    end

    % Calculate angles correspinding to coordinates of endpoints
    
    s = warning('off','MATLAB:divideByZero'); % suppress divide by 0 warnings
    endphi = atan(endpts(:,2)./endpts(:,1)) + pi*(endpts(:,1)<0);
    warning(s);
    endphi = unique(endphi);
    
    % One point is a nullset
    
    if(length(endphi)==1)
      continue
    end

    % Check whether or not any two successive phi form an interval inside
    % the triangle
    
    endphi = [endphi;endphi(1)+2*pi];
    phi0ind = [];
    
    inside = false;
    
    for j = 1:length(endphi)-1
      
      % no two consecutive intervals are inside
      
      if(inside)
        inside = true;
        continue
      end
      
      % Calculate the midpoint of the arc
      
      phi = 0.5*(endphi(j)+endphi(j+1));
      x = [cos(phi),sin(phi)];
      
      % x is inside the triangle iff it is in all of the
      % half-planes 'inside' all of the edges
      
      inside = true;
      
      for k = 1:3
        v1 = vcoord(mod(k-1,3)+1,:);
        v2 = vcoord(mod(k,3)+1,:);
        v3 = vcoord(mod(k+1,3)+1,:);
        
        w = v2-v1;
        n = [-w(2),w(1)];
        
        inside = inside && sign((v3-v1)*n')==sign((x-v1)*n');
      end
      
      if(inside)
        phi0ind = [phi0ind,j];
      end
       
    end % for loop (index j)
        
    % Do quadrature over each interval separately

    for j = 1:length(phi0ind)
      
      phi0 = endphi(phi0ind(j));
      phi1 = endphi(phi0ind(j)+1);
      
      if(phi0 == phi1)
        continue
      end
      
      % Transform quadrature points to [phi0,phi1]
      
      phi = phi0 + QuadRule.x(:)*(phi1-phi0);
      
      % Evaluate the function handle at phi
      
      FVal = FHandle(phi,Mesh.ElemFlag(i),varargin{:});
      
      % Transform quadrature points to triangle
      
      x_K = [cos(phi),sin(phi)];
      
      % Transform to reference triangle
      
      T = inv([vcoord(2,:)-vcoord(1,:);vcoord(3,:)-vcoord(1,:)]);
      
      x_Kref = (x_K-vcoord(ones(nPts,1),:))*T;
      
      % Evaluate shape functions
      
      N = shap_LFE(x_Kref);
      
      % Add contributions to global load vector
      
      L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1))*(phi1-phi0);
      L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2))*(phi1-phi0);
      L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3))*(phi1-phi0);
            
    end
    
  end % for loop over all elements
  
return