function I=integrate(n,m)

I=[];
for k=1:n
    for i=0:k
            j=k-i
            fhandle=@(x)x(:,1).^i.*x(:,2).^(j);

%             QuadRule=Duffy(TProd(gaulob(0,1,m)));
%             I_h=sum(QuadRule.w.*fhandle(QuadRule.x));

            QuadRule=NCC(m);
            I_h=[sum(QuadRule.w.*fhandle(QuadRule.x))];

           % QuadRule=Duffy(TProd(gauleg(0,1,10*m)));
            I_h=[I_h factorial(i)*factorial(j)/factorial(i+j+2)];
            I=[I; k i j I_h];        
    end
end

I