function UP_loc =STIMA_UPhp(Vertices, ElemFlag, EDofs, EDir, CDofs, EMass, VMass, QuadRule, VHandle, Shap,pmax); 
% STIMA_UPQFE  Element matrix for convection term using upwind
% quadrature rule 
%
%   LIELOC = STIMA_UPQFE(VERTICES,v, mass) computes the convection element 
%   matrix using bilinear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%    
%   MASS Mass of elements charing midpoints
%
%   Vhandle: handles velocity field
%    
%   Example:
%
%   lieloc = STIMA_UPQFE(Vertices,0,v);
%
%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
% Preallocate memory
  
  UP_loc = zeros(3+sum(EDofs)+CDofs,3+sum(EDofs)+CDofs);
  
 % initialize constant

Pts=size(QuadRule.w,1);

%UP_loc = zeros(6,6);
 
% Compute element mapping

a1 = Vertices(1,:);
a2 = Vertices(2,:);
a3 = Vertices(3,:);
  
bK = a1;
BK = [a2-bK; ...
        a3-bK];
det_BK = abs(det(BK));
inv_BK = inv(BK); 

% loop over all quadrature points
for i=1:Pts
    x=QuadRule.x(i,:);
    w=QuadRule.w(i);
    X=x*BK+a1;
    v=-VHandle(x*BK+a1);
    
    % compute element contribution
    loc_QR.x=x; loc_QR.w=w/det_BK;
    Shap_2D = shap_hp(loc_QR.x,pmax);
    Elem = STIMA_Conv_hp(Vertices, ElemFlag, EDofs, EDir, CDofs, loc_QR, Shap_2D, VHandle);
    % check location of quadrature point
    
    % first vertex
    if (x==[0 0])
        
        y = X + v;
        yhat = (y-bK)*inv_BK;
    
        z = sum(yhat);
        if(z > 0)
            xhat = yhat/z;
            if(xhat(1) >= 0 && xhat(2) >= 0)
                if (xhat(1)==0 || xhat(2)==0)   % two elements contribute to matrix
                        Elem=1/2 * Elem;
                end
                UP_loc=UP_loc+VMass(1)*Elem; 
            end
        end
    end
    
    
    % second vertex
    if  (x==[1 0])
        y = X + v;
        yhat = (y-bK)*inv_BK;
    
        z = 1-yhat(1);
        if(z > 0)
        xhat = [0 yhat(2)/z];
            if(xhat(2) >= 0 && xhat(2) <= 1)
                 if (xhat(2)==0 || xhat(1)== -xhat(2)) 
                     Elem=1/2 * Elem;   %two elements contribute to matrix
                 end  
                 UP_loc=UP_loc+VMass(2)*Elem; 
            end
        end    
    end
    
    
    % third vertex
    if (x==[0 1])
        y = X + v;
        yhat = (y-bK)*inv_BK;
  
        z = 1-yhat(2);
        if(z > 0)
            xhat = [yhat(1)/z 0];
            if(xhat(1) >= 0 && xhat(1) <= 1)
                if (xhat(1)==0 || xhat(1) == -xhat(2))
                  Elem=1/2 * Elem;   %two elements contribute to matrix 
                end
                UP_loc=UP_loc+VMass(3)*Elem; 
            end
        end      
    end
    
    % first edge (v1,v2)
    if (x(2)==0 && x(1)>0  && x(1)<1)
        yhat = (X+v-bK)*inv_BK;
        if(yhat(2) >= 0)
                if (yhat(2)==0)
                    Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+EMass(3)*Elem;
        end
    end
        
    % second edge (v2,v3)
    if (sum(x)==1 && x(1)>0  && x(1)<1)
        yhat = (X+v-bK)*inv_BK;
        if(sum(yhat) <= 1)
                if (sum(yhat)==1)
                     Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+EMass(1)*Elem;
        end
    end
    
    % third edge (v3,v1)
    if (x(1)==0 && x(2)>0  && x(2)<1)
        yhat = (X+v-bK)*inv_BK;
        if(yhat(1) >= 0)
                if (yhat(1)==0)
                     Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+EMass(2)*Elem;
        end
    end
    
end % loop over quadrature points

return