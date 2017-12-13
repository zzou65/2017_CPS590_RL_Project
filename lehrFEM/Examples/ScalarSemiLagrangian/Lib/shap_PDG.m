function shap = shap_PDG(x,p)
% SHAP_PDG Hierarchical shape functions.
%
%   SHAP = SHAP_PDG(X,P) computes the values and gradients of the
%   hierarchical shape functions on the reference element up to polynomial
%   degree P at the points X.
%
%   Example:
%
%   Shap = shap_PDG([0 0],10);

%   Copyright 2006-2011 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 s = shap_hp(x,p);
 
 shap = zeros(size(x,1),3+3*(p-1)+(p-2)*(p-1)/2);
 
 % vertex basis functions
 shap(:,1:3)=[s.vshap.shap{1} s.vshap.shap{2} s.vshap.shap{3}];
 
 % edge basis funtions
 offset = 3;
 for i = 0:p-2
     shap(:,offset + i*3 + 1) = s.eshap{1}.shap{i+1};
     shap(:,offset + i*3 + 2) = s.eshap{2}.shap{i+1};
     shap(:,offset + i*3 + 3) = s.eshap{3}.shap{i+1};
 end
 
 % cell basis functions
 offset = 3+3*(p-1);
 for i = 1:(p-2)*(p-1)/2
     shap(:,offset + i) = s.cshap.shap{i};
 end
 
 return