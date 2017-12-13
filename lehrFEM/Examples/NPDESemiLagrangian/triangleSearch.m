function [U,SearchComplete] = triangleSearch(X,P,UNode)
% triangleSearch supplied the Linear Interpolation value at a point if it
% present in the triangle formed by P or else supply 0
%
%   [U,SearchComplete] = triangleSearch(X,P,UNode) assembles stores Linear 
%   Interpolation value for point x in the triangle made by P and assign 1 
%   value to SearchComplete which is a search indicator or 0 if the point
%   lies outside search domain
%
%   by Rajdeep Deb
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%transformation of x from P triangle domain into [0 0;1 0;0 1] domain

P1 = X(1,:);
P2 = X(2,:);
P3 = X(3,:);
SearchComplete = 0;

bK = [P2 - P1 ; P3 - P1];

PT = ((P - P1)*inv(bK));

%Searching if x transformed exists in [0 0;1 0;0 1]

if PT(1,1) >= 0 && PT(1,2) >= 0 && sum(PT) <= 1
    N = shap_LFE(PT);
    U = UNode(1,1)*N(1,1) + UNode(2,1)*N(1,2) + UNode(3,1)*N(1,3);
    SearchComplete = 1;
else
    U = 0;
end;


