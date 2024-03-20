function ms=MS(pareto)
n=length(pareto);
f=zeros(n,2);

for i=1:n
    f(i,1)=pareto(i).Cost(1);
    f(i,2)=pareto(i).Cost(2);
end
minf1=min(f(:,1));
maxf1=max(f(:,1));
minf2=min(f(:,2));
maxf2=max(f(:,2));
ms=sqrt(((minf1-maxf1)^2)+((minf2-maxf2)^2));


end