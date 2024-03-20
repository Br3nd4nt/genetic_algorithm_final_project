function sns=SNS(MID,pareto)
n=length(pareto);
f=zeros(n,2);
M=0;
for i=1:n
    f(i,1)=pareto(i).Cost(1);
    f(i,2)=pareto(i).Cost(2);
end

for j=1:n
  ci=sqrt((f(j,1))^2+( f(j,2))^2);
 M=((MID-ci)^2)+M; 
end

sns=sqrt(M/(n-1));

end