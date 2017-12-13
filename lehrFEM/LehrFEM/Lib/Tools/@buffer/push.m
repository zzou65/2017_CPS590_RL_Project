function buf = push(buf,data)
% PUSH Push data inside FIFO buffer.
%
%   BUF = PUSH(BUF,DATA) push the matrix DATA inside the FIFO buffer BUF.
%
%   Example:
%
%   buf = open(buffer());
%   buf = push(buf,[1 2 3]);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  if(buf.open == buf.OPEN)
    try
    
      % Set position of file pointer to the end of the file
      
      if(buf.state == buf.READ_ONLY)
        fclose(buf.fid);
        buf.fid = fopen(buf.filename,'w');
        buf.state = buf.WRITE_ONLY;
      end
      fseek(buf.fid,0,'eof');
    
      % Write data field to buffer file
    
      [nrows,ncols] = size(data);
      fwrite(buf.fid,nrows,'double');
      fwrite(buf.fid,ncols,'double');
      if(isreal(data))
        fwrite(buf.fid,buf.REAL,'double');
        fwrite(buf.fid,data,'double');
        buf.nbytes = buf.nbytes + 8*(3+nrows*ncols);
      else
        fwrite(buf.fid,buf.COMPL,'double');
        fwrite(buf.fid,real(data),'double');
        fwrite(buf.fid,imag(data),'double');
        buf.nbytes = buf.nbytes + 8*(3+2*nrows*ncols);
      end
      buf.size = buf.size+1;
    catch
      error(ferror(fid));  
    end
  end
  
return