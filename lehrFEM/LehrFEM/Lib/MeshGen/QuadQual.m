function beta = QuadQual(Vertices)
% QUADQUAL Quality measure for quadrilateral elements.
%
%   BETA = QUADQUAL(VERTICES) computes the quality measure of a quadrilateral
%   element given by VERTICES.
%
%   Example:
%
%   Vertices = [0 0; 1 0; 1 1; 1 0];
%   beta = QuadQual(Vertices);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  AB = Vertices(2,:)-Vertices(1,:);
  AC = Vertices(3,:)-Vertices(1,:);
  BD_orth = (Vertices(4,:)-Vertices(2,:))*[0 -1; 1 0];
  
  beta = 0;
  if(abs(sum(AC.*BD_orth)) > 0) 
    QuadMid = Vertices(1,:) + sum(AB.*BD_orth)/sum(AC.*BD_orth)*AC;
    alpha = zeros(1,4);
    for i = 1:4
      AB = Vertices(i,:)-QuadMid;
      AC = Vertices(rem(i+4,4)+1,:)-QuadMid;
      alpha(i) = 2*sqrt(3)*(AB(1)*AC(2)-AB(2)*AC(1))/(sum(AB.^2)+sum(AC.^2)+sum((AC-AB).^2));
    end
    alpha = sort(alpha);
    beta = alpha(1)*alpha(2)/(alpha(3)*alpha(4));
  end
    
return