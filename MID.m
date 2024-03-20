function mid=MID(pareto)
Idealpoint=[0,1];
n=length(pareto);
f=zeros(n,2);
M=0;
for i=1:n
    f(i,1)=pareto(i).Cost(1);
    f(i,2)=pareto(i).Cost(2);
end

minf1=min(f(:,1));
maxf1=max(f(:,1));
minf2=min(f(:,2));
maxf2=max(f(:,2));

for j=1:n
   m1=(f(j,1)-Idealpoint(1))/(maxf1-minf1);
   m2=(f(j,2)-Idealpoint(2))/(maxf2-minf2); 
   M=sqrt(m1^2+m2^2)+M;
end
mid=M/n;

end