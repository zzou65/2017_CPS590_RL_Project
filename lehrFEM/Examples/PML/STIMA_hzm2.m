function Aloc = STIMA_hzm2(Vertices, flag, EHandle, QuadRule, varargin)
  
Aloc = zeros(3,3);

  l1x = Vertices(2,2)-Vertices(3,2); 
  l1y = Vertices(3,1)-Vertices(2,1);
  l2x = Vertices(3,2)-Vertices(1,2);
  l2y = Vertices(1,1)-Vertices(3,1); 
  l3x = Vertices(1,2)-Vertices(2,2);
  l3y = Vertices(2,1)-Vertices(1,1);


  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
   BK = [ P2 - P1 ; P3 - P1 ]; 
   det_BK = abs(det(BK));    
  nPoints = size(QuadRule.w,1);

  x = QuadRule.x*BK + ones(nPoints,1)*P1;
  
  Fval = EHandle(x,varargin{:});
  
  w = QuadRule.w;

  e = sum((Fval.*[w w w w]), 1);
  te(1,1) = e(1);
  te(1,2) = e(2);
  te(2,1) = e(3);
  te(2,2) = e(4);
  te = te/det_BK;

  Aloc(1,1) = (te*[l1x l1y]').'*[l1x l1y]';
  Aloc(1,2) = (te*[l1x l1y]').'*[l2x l2y]';
  Aloc(1,3) = (te*[l1x l1y]').'*[l3x l3y]';
  Aloc(2,2) = (te*[l2x l2y]').'*[l2x l2y]';
  Aloc(2,3) = (te*[l2x l2y]').'*[l3x l3y]';
  Aloc(3,3) = (te*[l3x l3y]').'*[l3x l3y]';
  Aloc(2,1) = (te*[l2x l2y]').'*[l1x l1y]';
  Aloc(3,1) = (te*[l3x l3y]').'*[l1x l1y]';
  Aloc(3,2) = (te*[l3x l3y]').'*[l2x l2y]';


return