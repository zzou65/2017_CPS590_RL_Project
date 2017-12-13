function g = gd_phi(x,varargin)
% implementation of the dirichlet boundary data for computing the cutoff
% function as the finite element solution of a differential equation

g = zeros(size(x,1),1);

outer = (x(:,1).^2 + x(:,2).^2)>0.95;
g(outer) = 1;