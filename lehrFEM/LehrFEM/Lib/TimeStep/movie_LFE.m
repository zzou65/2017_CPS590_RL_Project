function varargout = movie_LFE(buf,Mesh,nFrames,nRounds)
% MOVIE_LFE Generates a movie of a time dependent solution.
%
%   MOVIE_LFE(BUF,MESH,NFRAMES,NROUNDS) generates a movie of the time
%   dependent solution stored in the buffer variable BUF on the struct
%   MESH.
%   
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M
%                is equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh.
%      
%   F = MOVIE_LFE(BUF,MESH,NFRAMES,NROUNDS) also returns the structure of
%   the movie.
%
%   Example:
%
%   movie_LFE('sol_LFE.dat',Mesh,nSteps,100,5);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constant
  
  OFFSET = 0.05;
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Open up a tmeporary buffer
  
  tmp = open(buffer());
  
  % Compute color bar limits
  
  CMin = inf;
  CMax = -inf;
  nSteps = 0;
  while(~isempty(buf))
    [U,buf] = pop(buf);
    tmp = push(tmp,U);
    nSteps = nSteps+1;
    CMin = min(min(U),CMin);
    CMax = max(max(U),CMax);
  end
  if(CMin < CMax)          
    CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
  else
    CLim = [1-OFFSET 1+OFFSET]*CMin;   
  end
  nSteps = nSteps-1;
  
  % Record movie  
  
  step = 1;
  frame = 1;
  if(nFrames < nSteps+1)
    inc = floor((nSteps+1)/nFrames);
  else
    inc = 1;  
  end
  fig = figure('Name','Movie Recorder');
  while(~isempty(tmp))      
    [U,tmp] = pop(tmp);    
    if(mod(step,inc) == 0)
      view(2*[XMax YMax CMax]);
      axis ([XLim YLim CLim]);
      caxis([CMin CMax]);
      patch('Faces', Mesh.Elements, ...
            'Vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
            'CData', U, ...
            'Facecolor', 'interp', ...
            'Edgecolor', 'none');
      F(frame) = getframe;
      frame = frame+1;
      clf;
    end
    step = step+1;
  end
  close(fig);
  
  % Close the temporary buffer
  
  close(tmp);
  
  % Play the Movie
  
  fig = figure('Name','Movie Viewer');
  axis ([XLim YLim CLim]);
  movie(F,nRounds);
  if(nargout > 0)
    varargout{1} = F;  
  end
  
return