function eps=Epsilon(x,Sigma_Handle,varargin)
s=size(x,1);
eps=zeros(s,4);
for j=1:s
    kx=1;
    ky=1;
    if max(x(j,1))>0.75 
        kx=1-i*Sigma_Handle(x(j,1),0.75); end;
    if max(x(j,1))<-0.75 
        kx=1-i*Sigma_Handle(x(j,1),-0.75); end;
    if max(x(j,2))>0.75 
        ky=1-i*Sigma_Handle(x(j,2),0.75); end;
    if max(x(j,2))<-0.75 
        ky=1-i*Sigma_Handle(x(j,2),-0.75); end;
%     if max(x(j,1))>0.75 kx=1-i*50; end;
%     if max(x(j,1))<-0.75 kx=1-i*50; end;
%     if max(x(j,2))>0.75 ky=1-i*50; end;
%     if max(x(j,2))<-0.75 ky=1-i*50; end;
    eps(j,:)=[ky/kx 0 0 kx/ky];
end
return