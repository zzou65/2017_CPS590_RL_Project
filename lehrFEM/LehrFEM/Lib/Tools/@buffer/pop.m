function [data,buf] = pop(buf);
% POP Pop data from FIFO buffer file.
%
%   [DATA,BUF] = POP(BUF,POS) pops the data matrix DATA from the FIFO
%   buffer BUF
%
%   Example:
%
%   buf = open(buffer());
%   buf = push(buf,[1 2 3]);
%   [data,buf] = pop(buf);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  if(buf.open == buf.OPEN && ~isempty(buf))

    % Set position of file pointer
  
    try
      if(buf.state == buf.WRITE_ONLY)
        fclose(buf.fid);
        fid = fopen(buf.filename,'r');
        buf.state = buf.READ_ONLY;
      end
      fseek(buf.fid,buf.pos,BOF);
    catch
      error(ferror(buf.fid));  
    end
  
    % Extract real data from file
  
    try
      nrows = fread(buf.fid,1,'*double');
      ncols = fread(buf.fid,1,'*double');
      type = fread(buf.fid,1,'*double');
      if(type == buf.REAL)
        data = fread(buf.fid,nrows*ncols,'*double');
        buf.pos = buf.pos + 8*(3+nrows*ncols);
      else
        data = fread(buf.fid,nrows*ncols,'*double') + ...
               i*fread(buf.fid,nrows*ncols,'*double');
        buf.pos = buf.pos + 8*(3+2*nrows*ncols);
      end
    catch
      error(ferror(buf.fid));  
    end
    buf.size = buf.size-1;
    data = reshape(data,nrows,ncols);
  
  else  
    data = [];
  end
  
return