function dist = dist_diff(dist_1,dist_2)
% DIST_DIFF  Signed distance for the difference between two domains.
%
%   DIST = DIST_DIFF(DIST_1,DIST_2) computes the value of the signed distance
%   function for the difference between two domains with signed distance
%   functions DIST_1 and DIST_2. 
%
%   Example:
%
%   dist = dist_diff(dist_1,dist_2);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  dist = max(dist_1,-dist_2);
  
return
