function dist = dist_isect(dist_1,dist_2)
% DIST_ISECT Signed distance for the intersection of two domains.
%
%   DIST = DIST_ISECT(DIST_1,DIST_2) computes the value of the signed distance
%   function for the intersection of two domains with signed distance functions
%   DIST_1 and DIST_2. 
%
%   Example:
%
%   dist = dist_isect(dist_1,dist_2);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  dist = max(dist_1,dist_2);
  
return