function v = vel_circ(x,varargin)
%VEL_CIRC circular velocity
%   Detailed explanation goes here

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

m = [0.5,0.5];
x0 = x-m(ones(size(x,1),1),:);
v = [-x0(:,2),x0(:,1)];