function [] = display(buf)
% DISPLAY Display a FIFO buffer.
%
%   DISPLAY(BUF) displays a FIFO buffer BUF on the screen.
% 
%   Example:
%
%   display(buf)

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Display buffer on screen

  if(buf.open == buf.OPEN)
    fprintf('\n');
    fprintf([' ' inputname(1) ' =\n']);
    fprintf('\n');
    fprintf(['   buffer [open, ' int2str(size(buf)) ' elements]\n']);
    fprintf('\n');
  else
    fprintf('\n');
    fprintf([' ' inputname(1) ' =\n']);
    fprintf('\n');
    fprintf(['   buffer [closed]\n']);
    fprintf('\n');
  end
  
return