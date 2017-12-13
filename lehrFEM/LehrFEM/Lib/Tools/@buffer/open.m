function buf = open(buf,varargin)
% OPEN Open FIFO buffer.
%
%   BUF = OPEN(BUF) open up a FIFO buffer BUF.
%
%   BUF = OPEN(BUF,PATH) open up a FIFO buffer in the directoy specified by
%   PATH. 
%
%   Example:
%
%   buf = open(buffer());
  
%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize buffer counter

  persistent count;
  if(isempty(count))
    count = 0;  
  end
  
  if(buf.open == buf.CLOSED)
  
    % Generate filename for current buffer
    
    if(nargin > 1)
      filename = [varargin{1} '/buffer_' int2str(count) '.dat'];
    else
      filename = ['buffer_' int2str(count) '.dat']; 
    end
    
    % Create empty buffer file
  
    try
      fid = fopen(filename,'w');
      fclose(fid);
    catch
      error(ferror(fid)); 
    end

    % Open buffer file for read and write access
  
    buf.open = buf.OPEN;
    buf.filename = filename;
    try  
      buf.fid = fopen(buf.filename,'a');
    catch
      error(ferror(fid));  
    end
    buf.state = buf.WRITE_ONLY;
    buf.pos = 0;
    buf.size = 0;
    buf.nbytes = 0;
    
    % Update open buffer counter
    
    count = count+1;
  end
  
return