function [] = save_Mesh(Mesh,CoordFile,ElemFile)
% SAVE_MESH Save mesh to file.
%
%   MESH = SAVE_MESH(MESH,COORDFILE,ELEMFILE) saves a mesh to the files
%   COORDFILE (list of vertices) and ELEMFILE (list of elements).
%
%   The struct MESH must at least contain the followng fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   save_Mesh(Mesh,'Coordinates.dat','Elements.dat');

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize global constants
  
  global C I_PREC D_PREC DELIM NL;
  
  C = '%%';           % Character specifying Matlab comments
  I_PREC = '%-5d';    % Field width for integers
  D_PREC = '%+2.6e';  % Field width for doubles
  DELIM = '  ';       % Field delimiter
  NL = '\n';          % Newline character
  
  % Save mesh data to files
  
  save_Coordinates(Mesh.Coordinates,CoordFile);
  save_Elements(Mesh.Elements,ElemFile);
  
return

%%% Saves coordinates to files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = save_Coordinates(Coordinates,File)
% SAVE_COORDINATES Save vertex coordinates to file.
%
%   SAVE_COORDINATES(COORDINATES,FILE) saves the vertex coordinates to the
%   .dat file FILE.
%
%   Example:
%
%   save_Coordinates(Coordinates,'Coordinates.dat');
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  global I_PREC D_PREC DELIM NL;

  % Open file with writing permissions

  fid = fopen(File,'w');
  if(fid < 2)
    error(['Could not open file ',File]);  
  end

  % Create file header

  Str = Header(1);
  fprintf(fid,Str);
  fprintf(fid,'\n');
  
  % Write data to file
  
  nCoordinates = size(Coordinates,1);
  for i = 1:nCoordinates
    fprintf(fid,[sprintf(I_PREC,i), DELIM, ...
                 sprintf(D_PREC,Coordinates(i,1)), DELIM, ...
                 sprintf(D_PREC,Coordinates(i,2)), NL]);
  end
   
  
  % Close file
  
  fclose(fid);
  
return

%%% Saves elements to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = save_Elements(Elements,File)
% SAVE_ELEMENTS Save elements to file.
%
%   SAVE_ELEMENTS(ELEMENTS,FILE) saves the elements of a mesh to the .dat
%   file FILE.
%
%   Example:
%
%   save_Elements(Elements,'Elements.dat');
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  global I_PREC D_PREC DELIM NL;

  % Open file with writing permissions

  fid = fopen(File,'w');
  if(fid < 2)
    error(['Could not open file ',File]);  
  end

  % Create file header

  Str = Header(2);
  fprintf(fid,Str);
  fprintf(fid,'\n');
  
  % Write data to file
  
  [nElements,nVert] = size(Elements);
  if(nVert == 3)
    
    % Save triangular elements  
      
    for i = 1:nElements
      fprintf(fid,[sprintf(I_PREC,i), DELIM, ...
                   sprintf(I_PREC,Elements(i,1)), DELIM, ...
                   sprintf(I_PREC,Elements(i,2)), DELIM, ...
                   sprintf(I_PREC,Elements(i,3)), NL]);  
    end 
  else
      
    % Save quadrilateral elements  
      
    for i = 1:nElements
      fprintf(fid,[sprintf(I_PREC,i), DELIM, ...
                   sprintf(I_PREC,Elements(i,1)), DELIM, ...
                   sprintf(I_PREC,Elements(i,2)), DELIM, ...
                   sprintf(I_PREC,Elements(i,3)), DELIM, ...
                   sprintf(I_PREC,Elements(i,4)), NL]);  
    end  
  end
    
  % Close file
  
  fclose(fid);
  
return

%%% Generates header of mesh files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Str = Header(Type)
% HEADER Generates header strings.
%
%   STR = HEADER(TYPE) generates header strings according to the
%   integer TYPE:
%    1 Header string for coordinate files
%    2 Header string for element files
%
%   Example:
%
%   Str = Header(1);
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  global C NL;

  switch(Type)
    case 1  
      Str = [C, ' List of vertices', NL, ...
             C, NL, ...
             C, '   User  :  ', upper(getenv('USER')), NL, ...
             C, '   Date  :  ', datestr(now,0), NL, ...
             C, NL, ...
             C, '   The first column contains the global index of each vertex.', NL, ...
             C, NL];
    case 2
      Str = [C, ' List of elements', NL, ...
             C, NL, ...
             C, '   User  :  ', upper(getenv('USER')), NL, ...
             C, '   Date  :  ', datestr(now,0), NL, ...
             C, NL, ...
             C, '   The first column contains the global index of each element.', NL, ...
             C, NL];
  end
  
return
