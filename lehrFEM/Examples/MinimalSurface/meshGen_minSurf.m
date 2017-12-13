% meshGen_minSurf
%   Run script for creating graded square meshes for the minimal surface
%   problem. See GRADED_SQUARE for more info

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc

% Initialize constants

Beta = 2.7;     % Grading factor
N = 4;          % Number of meshes to create
nNodes = 5;
h = 0.06;

for i = 1:N,
    Mesh = Sqr_graded(h);
    save_Mesh(Mesh, ...
              ['Coord_Sqr_gr_' int2str(i) '.dat'], ...
              ['Elem_Sqr_gr_' int2str(i) '.dat']);
    close all
    h = h/2;
end