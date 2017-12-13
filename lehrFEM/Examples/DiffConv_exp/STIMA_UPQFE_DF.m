function UP_loc = STIMA_UPQFE_DF(Vertices, Flag, vHandle, varargin)
% STIMA_UPQFE  Element matrix for convection term using upwind
% quadrature rule with weights at vertices, midpoints and barycenter
%
%   LIELOC = STIMA_UPQFE(VERTICES,v) computes the convection element 
%   matrix using bilinear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Vhandle: handles velocity field
%    
%
%   Example:
%
%   lieloc = STIMA_UPQFE(Vertices,0,v);
%
%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% initialize constant

UP_loc = zeros(6,6);
 
% Compute element mapping

a1 = Vertices(1,:);
a2 = Vertices(2,:);
a3 = Vertices(3,:);
  
bK = a1;
BK = [a2-bK; ...
        a3-bK];
inv_BK = inv(BK); 
  
% Compute gradient of local shape functions
QuadRule.x=[0 0; 1 0; 0 1; 0.5 0.5; 0 0.5; 0.5 0]; 

Pts=size(QuadRule.x,1);

N=shap_QFE(QuadRule.x);  
gradN=grad_shap_QFE(QuadRule.x);
for i=1:2:12
      gradN(:,[i i+1])=gradN(:,[i i+1])*inv_BK';
end

% loop over all quadrature points
for i=1:Pts
    x=QuadRule.x(i,:);
    X=x*BK+a1;
    v=-vHandle(x*BK+a1);
    
    % check location of quadrature point
    
    % first vertex
    if (x==[0 0])
        
        y = X + v;
        yhat = (y-bK)*inv_BK;
    
        z = sum(yhat);
        if(z > 0)
            xhat = yhat/z;
            if(xhat(1) >= 0 && xhat(2) >= 0)
                Elem =N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                if (xhat(1)==0 || xhat(2)==0)   % two elements contribute to matrix
                        Elem=1/2 * Elem;
                end
                UP_loc=UP_loc+Elem; 
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
                Elem = N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                 if (xhat(2)==0 || xhat(1)== -xhat(2)) 
                     Elem=1/2 * Elem;   %two elements contribute to matrix
                 end  
                 UP_loc=UP_loc+Elem; 
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
                Elem =N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                if (xhat(1)==0 || xhat(1) == -xhat(2))
                  Elem=1/2 * Elem;   %two elements contribute to matrix 
                end
                UP_loc=UP_loc+Elem; 
            end
        end      
    end
    
    % first edge
    if (x(2)==0 && x(1)>0  && x(1)<1)
        yhat = (X+v-bK)*inv_BK;
        if(yhat(2) >= 0)               
                Elem = N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                if (yhat(2)==0)
                    Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+Elem;
        end
    end
        
    % second edge
    if (sum(x)==1 && x(1)>0  && x(1)<1)
        yhat = (X+v-bK)*inv_BK;
        if(sum(yhat) <= 1)
                Elem = N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                if (sum(yhat)==1)
                     Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+Elem;
        end
    end
    
    % third edge
    if (x(1)==0 && x(2)>0  && x(2)<1)
        yhat = (X+v-bK)*inv_BK;
        if(yhat(1) >= 0)
               Elem = N(i,:)'*[-gradN(i,[1 2])*v', -gradN(i,[3 4])*v', ....
                            -gradN(i,[5 6])*v', -gradN(i,[7 8])*v', ....
                            -gradN(i,[9 10])*v', -gradN(i,[11 12])*v'];
                if (yhat(1)==0)
                     Elem=Elem/2;    %two elements contribute to matrix 
                end
                UP_loc=UP_loc+Elem;
        end
    end
    
end % loop over quadrature points

return
