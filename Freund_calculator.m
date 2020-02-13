%Recall: MATLAB does array coordinates (x,y) --> (column, row)
%           rather than standard (row, column)

n = 4;
k = 30;

cRes = 100;
c = linspace(0.1,2,cRes);

pExact = zeros(n,n,cRes);
sigma = zeros(n);

for p=1:n
    for q=1:n
        sigma(p,q) = n/(2*pi^2*sqrt(k)*(p^2+q^2));
        for r=1:cRes
            pExact(p,q,r) = sqrt(2/pi)* exp(-c(r)^2/(2*sigma(p,q)^2)) / ...
                sigma(p,q)*erf(c(r)/(sqrt(2)*sigma(p,q)));
        end
    end
end

fExact = squeeze(sum(pExact,[1 2]));
fAppx = n^2./c;

figure(1)
hold on
plot(c,fAppx,'-r')
plot(c,fExact,'-b')
xlabel('F')
ylabel('c')
hold off