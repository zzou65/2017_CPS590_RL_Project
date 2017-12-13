function d=getData_EvHs();

% Meshdata
data.boundtype=[-1 -1 -2 -2];
data.Coordinates=[0 0; 1 0; 1 1; 0 1];
data.Elements=[1 2 4;2 3 4];

EX1=@(x,a) exp(2/a*(x-1));
DEX1=@(x,a) 2/a*exp(2/a*(x-1));
DDEX1=@(x,a) 2/a*2/a*exp(2/a*(x-1));
DDDEX1=@(x,a) 2/a*2/a*2/a*exp(2/a*(x-1));
EX2=@(x,a) exp(3/a*(x-1));
DEX2=@(x,a) 3/a*exp(3/a*(x-1));
DDEX2=@(x,a) 3/a*3/a*exp(3/a*(x-1));
DDDEX2=@(x,a) 3/a*3/a*3/a*exp(3/a*(x-1));

Hs=@(x,a)...
    x(:,1).*x(:,2).^2-...
    x(:,2).^2.*EX1(x(:,1),a)-...
    x(:,1).*EX2(x(:,2),a)+...
    EX1(x(:,1),a).*EX2(x(:,2),a);
D1Hs=@(x,a)...
    x(:,2).^2-...
    x(:,2).^2.*DEX1(x(:,1),a)-...
    EX2(x(:,2),a)+...
    DEX1(x(:,1),a).*EX2(x(:,2),a);
D2Hs=@(x,a)...
    2*x(:,1).*x(:,2)-...
    2*x(:,2).*EX1(x(:,1),a)-...
    x(:,1).*DEX2(x(:,2),a)+...
    EX1(x(:,1),a).*DEX2(x(:,2),a);
D11Hs=@(x,a)...
    -x(:,2).^2.*DDEX1(x(:,1),a)+...
    DDEX1(x(:,1),a).*EX2(x(:,2),a);
D12Hs=@(x,a)...
    2*x(:,2)-...
    2*x(:,2).*DEX1(x(:,1),a)-...
    DEX2(x(:,2),a)+...
    DEX1(x(:,1),a).*DEX2(x(:,2),a);
D22Hs=@(x,a)...
    2*x(:,1)-...
    2*EX1(x(:,1),a)-...
    x(:,1).*DDEX2(x(:,2),a)+...
    EX1(x(:,1),a).*DDEX2(x(:,2),a);
D111Hs=@(x,a)...
    -x(:,2).^2.*DDDEX1(x(:,1),a)+...
    DDDEX1(x(:,1),a).*EX2(x(:,2),a);
D112Hs=@(x,a)...
    -2*x(:,2).*DDEX1(x(:,1),a)+...
    DDEX1(x(:,1),a).*DEX2(x(:,2),a);
D122Hs=@(x,a)...
    2*ones(size(x,1),1)-...
    2*DEX1(x(:,1),a)-...
    DDEX2(x(:,2),a)+...
    DEX1(x(:,1),a).*DDEX2(x(:,2),a);
D222Hs=@(x,a)...
    -x(:,1).*DDDEX2(x(:,2),a)+...
    EX1(x(:,1),a).*DDDEX2(x(:,2),a);

Ev1=@(x,a) -D2Hs(x,a);
D1Ev1=@(x,a) -D12Hs(x,a);
D2Ev1=@(x,a) -D22Hs(x,a); 
D11Ev1=@(x,a) -D112Hs(x,a);
D12Ev1=@(x,a) -D122Hs(x,a); 
D22Ev1=@(x,a) -D222Hs(x,a); 

Ev2=@(x,a) D1Hs(x,a);
D1Ev2=@(x,a) D11Hs(x,a);
D2Ev2=@(x,a) D12Hs(x,a);
D11Ev2=@(x,a) D111Hs(x,a);
D12Ev2=@(x,a) D112Hs(x,a);
D22Ev2=@(x,a) D122Hs(x,a);
 
% Ev1=@(x,a)sin(pi.*x(:,1)).*(1-x(:,2));
% D1Ev1=@(x,a)pi*(1-x(:,2)).*cos(pi.*x(:,1));
% D2Ev1=@(x,a)-sin(pi.*x(:,1));
% D11Ev1=@(x,a)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
% D12Ev1=@(x,a)-pi.*cos(pi.*x(:,1));
% D22Ev1=@(x,a)0*x(:,1);
% 
% Ev2=@(x,a) (1-x(:,2).^2).*(1-x(:,1).^2);
% D1Ev2=@(x,a)(1-x(:,2).^2).*(-2*x(:,1));
% D2Ev2=@(x,a)(-2*x(:,2)).*(1-x(:,1).^2);
% D11Ev2=@(x,a)(1-x(:,2).^2).*(-2);
% D12Ev2=@(x,a)(-2*x(:,2)).*(-2*x(:,1));
% D22Ev2=@(x,a)(-2).*(1-x(:,1).^2);

%V1=@(x)2*ones(size(x,1),1);
%V2=@(x)3*ones(size(x,1),1);
 v1=2;
 v2=3;
 V1=@(x)v1*x(:,2).*x(:,1);
 V2=@(x)v2*x(:,2).*x(:,2);

Curlfv = @(x,a) a*D11Hs(x,a)+a*D22Hs(x,a)+V1(x)*D1Hs(x,a)+V2(x)*D2Hs(x,a)-Hs(x,a);
%fv=@(x,a) [D2Hs(x,a)+V2(x)*Hs(x,a) -D1Hs(x,a)-V2(x)*Hs(x,a)];
fv=@(x,a) [-a*D22Ev1(x,a)+a*D12Ev2(x,a)-V2(x).*(-D2Ev1(x,a)+D1Ev2(x,a))+Ev1(x,a) ...
                 a*D12Ev1(x,a)-a*D11Ev2(x,a)-V1(x).*(D2Ev1(x,a)-D1Ev2(x,a))+Ev2(x,a)];
% velocity field
data.V_Handle=@(x,varargin)[V1(x) V2(x)];
data.JacV_Handle=@(x,t)zeros(size(x,1),4);
% Hs
data.Hs_Handle=@(x,d,a)Hs(x,a);
% grad Hs
data.GRAD_Hs_Handle=@(x,d,a)[D1Hs(x,a) D2Hs(x,a)];
% curl fv
data.CurlFv_Handle=@(x,d,a)Curlfv(x,a);
% Ev
data.Ev_Handle=@(x,d,a)a*[Ev1(x,a) Ev2(x,a)];
% -curl Ev
data.CURL_Ev_Handle=@(x,d,a)-a*D1Ev2(x,a)+a*D2Ev1(x,a);
% - v x curl Ev
data.VCURL_Ev_Handle=@(x,d,a)-a*[(D2Ev1(x,a)-D1Ev2(x,a)).*V2(x) (-D2Ev1(x,a)+D1Ev2(x,a)).*V1(x)];
% v  Ev
data.V_Ev_Handle=@(x,d,a)-a*(V1(x)*Ev1(x,a)+V2(x)*Ev2(x,a));
% fv
data.Fv_Handle=@(x,d,a)a*fv(x,a);

d=data;

return