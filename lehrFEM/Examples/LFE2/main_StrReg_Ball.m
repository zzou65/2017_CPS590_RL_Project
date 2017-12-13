% Run script for strong regularization LFE2 finite element solver on the 
% disk domian.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

  DHANDLE = @dist_circ;  % Signed distance function for a circle
  Center = [0 0];             % Center of circle
  R = 1;                 % Radius of circle
  NREFS = 2;             % Total number mesh refinements
  U_Handle = @(x,varargin)ones(size(x,1),1);
  F_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))]+[sin(pi*x(:,2)) sin(pi*x(:,1))];
  GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
%   F_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];
%   GD_Handle = @(x,varargin)[-x(:,2) x(:,1)];

  % Initialize mesh

  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  BdFlags = -1;
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = BdFlags;
  Mesh = refine_REG(Mesh,DHANDLE,Center,R);
  
  % Do initial mesh refinements

  for i = 1:NREFS
      
      Mesh = refine_REG(Mesh,DHANDLE,Center,R);

      % Assemble Curl-curl matrix, MASS matrix and load vector

      nCoordinates = size(Mesh.Coordinates,1);
      [IC,JC,C] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
      [ID,JD,D] = assemMat_LFE2(Mesh,@STIMA_Div_LFE2,U_Handle,P7O6());
      [U,g,FreeDofs,IB,JB,B] = assemDir_StrRegLFE2(Mesh,BdFlags,GD_Handle);
      A = sparse([IC;ID;IB+2*nCoordinates;JB],[JC;JD;JB;IB+2*nCoordinates],[C;D;B;B]);
      l = assemLoad_LFE2(Mesh,P7O6(),F_Handle);
      L = [l;g];
      L = L - A*U;
      
      % Solve the system

      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
      u = U(1:2*nCoordinates);

      % Plot the solution

      plot_DomBd(Mesh,'as')
      hold on
      plot_LFE2(u,Mesh);

  end

  % Clear memory

  clear all
  
  
  