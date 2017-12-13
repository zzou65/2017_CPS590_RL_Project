function eta = circref(U,Mesh)
%CIRCREF refinement towards S1

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  numelem = size(Mesh.Elements,1);

  eta = zeros(numelem,1);
  
  for i = 1:numelem
  
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    vcoord = Mesh.Coordinates(vidx,:);
    
    % Find points on intersection of unit circle and boundary of element
    
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
          eta(i) = 1;
          break
        end 
      else
        discr = b^2-4*a*c;

        if(discr == 0)
          t = -b/(2*a);
          if(t>=0 && t<=1)
            eta(i) = 1;
            break
          end
        elseif(discr > 0)
          t = (-b-sqrt(discr))/(2*a);
          if(t>=0 && t<=1)
            eta(i) = 1;
            break
          end
          t = (-b+sqrt(discr))/(2*a);
          if(t>=0 && t<=1)
            eta(i) = 1;
            break
          end
        end
      end
      
    end % for loop over edges of triangle (index j)
    
  end

return