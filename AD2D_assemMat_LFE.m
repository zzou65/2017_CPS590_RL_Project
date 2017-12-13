function varargout = AD2D_assemMat_LFE(Mesh,EHandle,EleleMent_Props,varargin)

  nElements = size(Mesh.Elements,1);
  % Preallocate memory
  I = zeros(9*nElements,1);
  J = zeros(9*nElements,1);
  A = zeros(9*nElements,1);
  
  % Check for element flags
  if(isfield(Mesh,'ElemFlag')),
    flags = Mesh.ElemFlag; 
  else
    flags = zeros(nElements,1);
  end

  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    EleProp = EleleMent_Props(i);
    % Extract vertices of current element
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
      
    % Compute element contributions
    Aloc = EHandle(Vertices,flags(i),varargin{:});
   
    % mutiply by EleProp
    Aloc = EleProp * Aloc;
    
    % Add contributions to stiffness matrix
    I(loc) = set_Rows(idx,3);
    J(loc) = set_Cols(idx,3);
    A(loc) = Aloc(:);
    loc = loc+9;
  end
  
  
  % Assign output arguments
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);      
  end
  
return