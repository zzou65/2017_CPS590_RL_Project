function Mesh = load_Mesh(CoordFile,ElemFile)
% LOAD_MESH Load mesh from file.
%
%   MESH = LOAD_MESH(COORDFILE,ELEMFILE) loads a mesh from the files COORDFILE
%   (list of vertices) and ELEMFILE (list of elements).
%
%   The struct MESH contains the followng fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   Mesh = load_Mesh('Coordinates.dat','Elements.dat');

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Load mesh from files

  Mesh.Coordinates = load_Coordinates(CoordFile);
  Mesh.Elements = load_Elements(ElemFile);
  
return

%%%

function Coordinates = load_Coordinates(File)
% LOAD_COORDINATES Load vertex coordinates from file.
%
%   COORDINATES = LOAD_COORDINATES(FILE) loads the vertex coordinates from
%   the .dat file FILE.
%
%   Example:
%
%   Coordinates = load_Coordinates('Coordinates.dat');
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  Coordinates = load(File);
  Coordinates(:,1) = [];

return

%%%

function Elements = load_Elements(File)
% LOAD_ELEMENTS Load elements from a file.
%
%   ELEMENETS = LOAD_ELEMENTS(FILE) load the elements of a mesh from the
%   .dat file FILE.
%
%   Example:
%
%   Elements = load_Elements('Elements.dat');
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  Elements = load(File);
  Elements(:,1) = [];
    
return
