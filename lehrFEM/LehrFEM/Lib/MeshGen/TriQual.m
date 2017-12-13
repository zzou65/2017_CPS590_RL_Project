function alpha = TriQual(Vertices)
% TRIQUAL Quality measure for triangular elements.
%
%   ALPHA = TRIQUAL(VERTICES) computes the quality measure of a triangular
%   element given by VERTICES.
%
%   Example:
%
%   Vertices = [0 0; 1 0; 0 1];
%   alpha = TriQual(Vertices);
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  AB = Vertices(2,:)-Vertices(1,:);
  AC = Vertices(3,:)-Vertices(1,:);
  alpha = 2*sqrt(3)*(AB(1)*AC(2)-AB(2)*AC(1))/(sum(AB.^2)+sum(AC.^2)+sum((AC-AB).^2));
  
return