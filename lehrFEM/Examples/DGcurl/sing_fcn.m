function fval = sing_fcn(x)
fval = 0*ones(size(x,1),1);
for j = 1:size(x,1)
    if norm(x(j,:))>0
        fval(j) = 2/3*(x(j,1).^2+x(j,2).^2).^(-1/6);
    end
end