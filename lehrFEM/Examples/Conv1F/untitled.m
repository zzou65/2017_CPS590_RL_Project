m=min(real(eig(J([1,0]))));
for i=-500:500
    for j=-500:500
        m=min(m,min(real(eig(J([i/500,j/500])))));
    end
    m
end
        