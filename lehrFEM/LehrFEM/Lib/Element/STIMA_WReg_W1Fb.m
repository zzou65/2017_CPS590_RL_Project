function Aloc = STIMA_WReg_W1F(Vertices,ElemInfo,QuadRule,varargin)
% Using QuadRule for the future work(space dependent version)

    % Initialize constant
    
    nGuass = size(QuadRule.w,1);
    
    % Preallocate memory
    
    Aloc = zeros(3,3);
    N_W1F = shap_W1F(QuadRule.x);
    grad_N = grad_shap_LFE(QuadRule.x);
    
    % Compute element mapping

    P1 = Vertices(1,:);
    P2 = Vertices(2,:);
    P3 = Vertices(3,:);
    bK = P1;
    BK = [P2-bK;P3-bK];
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    TK = transpose(inv_BK);
    
    % Compute element entry
    
    N(:,1:2) = N_W1F(:,1:2)*TK;
    N(:,3:4) = N_W1F(:,3:4)*TK;
    N(:,5:6) = N_W1F(:,5:6)*TK;
    G(:,1:2) = grad_N(:,1:2)*TK;
    G(:,3:4) = grad_N(:,3:4)*TK;
    G(:,5:6) = grad_N(:,5:6)*TK;
    
    Aloc(1,1) = sum(QuadRule.w.*sum(N(:,1:2).*G(:,1:2),2))*det_BK;
    Aloc(1,2) = sum(QuadRule.w.*sum(N(:,1:2).*G(:,3:4),2))*det_BK;
    Aloc(1,3) = sum(QuadRule.w.*sum(N(:,1:2).*G(:,5:6),2))*det_BK;
    Aloc(2,1) = sum(QuadRule.w.*sum(N(:,3:4).*G(:,1:2),2))*det_BK;
    Aloc(2,2) = sum(QuadRule.w.*sum(N(:,3:4).*G(:,3:4),2))*det_BK;
    Aloc(2,3) = sum(QuadRule.w.*sum(N(:,3:4).*G(:,5:6),2))*det_BK;
    Aloc(3,1) = sum(QuadRule.w.*sum(N(:,5:6).*G(:,1:2),2))*det_BK;
    Aloc(3,2) = sum(QuadRule.w.*sum(N(:,5:6).*G(:,3:4),2))*det_BK;
    Aloc(3,3) = sum(QuadRule.w.*sum(N(:,5:6).*G(:,5:6),2))*det_BK;
    
return