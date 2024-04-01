function model = model()
%% selection of the size of test problems:
%%%%%%%%%% P1

m=23; % Number of suppliers 
p=13; % Number of potential plants 

n=16; % Number of potential distributers 
i=21; % Number of customers 
j=16;  % Number of potential recovering centers 
d=11;  % Number of potential remanufacturing ecnters 
r=10;  % Number of potential recycling centers 
s=6;  % Number of disposal centers 

% m=3; % Number of suppliers 
% p=6; % Number of potential plants 
% n=12; % Number of potential distributers 
% i=20; % Number of customers 
% d=4;  % Number of potential remanufacturing ecnters 
% r=3;  % Number of potential recycling centers 
% j=3;  % Number of potential recovering centers 
% s=4;  % Number of disposal centers 



c=10;  % Number of scenarios 

model.m=m; 
model.p=p; 
model.n=n; 
model.i=i; 
model.d=d; 
model.r=r; 
model.j=j; 
model.s=s; 
model.c=c; 

%% parameters of CLSC Problem under uncertainty 
pbc=zeros(1,c); 
for cc=1:c
    pbc(cc)=(1/c); 
end 
model.pbc=pbc; % probability of scenarios 

accni=randi([0, 12],n,i,c)+rand; 
model.accni=accni; % assignment cost 

mcpc=randi([5, 32],p,c)+rand; 
model.mcpc=mcpc;  % varibale cost of manufacturing 

pcmc=randi([5, 32],p,c)+rand; 
model.pcmc=pcmc;  % Price of manufactured products 

% Transportation costs
tcmp=randi([3, 8],m,p); 
model.tcmp=tcmp; 

tcpn=randi([3, 8], p, n); 
model.tcpn=tcpn; 

tcsm=randi([3, 8], s, m); 
model.tcsm=tcsm; 

tcjs=randi([3, 8],j, s); 
model.tcjs=tcjs; 

tcjd=randi([3, 8], j, d); 
model.tcjd=tcjd; 

tcdp=randi([3, 8], d, p);
model.tcdp=tcdp; 

tcjr=randi([3, 8], j, r);
model.tcjr=tcjr; 

tcrn=randi([3, 8], r, n);
model.tcrn=tcrn; 

tcij=randi([3, 8], i, j); 
model.tcij=tcij; 

% Maximum desired number of facilities
maxd=randi([round(d/2)-1, round(d/2)]); 
model.maxd=maxd; 

maxm=randi([round(m/2)-1, round(m/2)]); 
model.maxm=maxm; 

maxn=randi([round(n/2)-1, round(n/2)]); 
model.maxn=maxn; 

maxp=randi([round(p/2)-1, round(p/2)]); 
model.maxp=maxp; 

maxr=randi([round(r/2)-1, round(r/2)]); 
model.maxr=maxr; 

maxs=randi([round(s/2)-1, round(s/2)]); 
model.maxs=maxs; 

maxj=randi([round(j/2)-1, round(j/2)]); 
model.maxj=maxj; 

dci=normrnd(randi([100, 260]), 20, c, i); % Demand of customers 
model.dci=dci;

% Fraction rates 
bp=rand(1, p)*0.03; 
model.bp=bp; 
alphai=0.3+rand(1,i)*0.4; 
model.alphai=alphai; 

% Profits from recovery activities
scj=randi([5, 10], 1, j);
model.scj=scj; 
scd=randi([5, 10], 1, d);
model.scd=scd; 
scr=randi([5, 10], 1, r); 
model.scr=scr; 

% Rates for retunred products 
a1=rand*0.4; 
a2=rand*0.4; 
a3=1-a1-a2; 
model.a1=a1;
model.a2=a2;
model.a3=a3; 

% Capacity of facilities 
landa=0; 
for ii=1:i
    for cc=1:c
        landa=landa+dci(cc,ii); 
    end
end
landa=landa/c; 
betta=rand; 
gamma=rand; 
pn=randi([round((1-gamma)*((landa*(1-betta))/maxn)), round((1+gamma)*((landa*(1+betta))/maxn))],1,n); 
model.pn=pn; 
landa=0; 
for ii=1:i
    for cc=1:c
        landa=landa+dci(cc,ii)*alphai(ii); 
    end
end
pj=randi([round((1-gamma)*((landa*(1-betta))/maxj)), round((1+gamma)*((landa*(1+betta))/maxj))],1,j);
model.pj=pj; 
pd=randi([round((1-gamma)*((a1*landa*(1-betta))/maxd)), round((1+gamma)*((a1*landa*(1+betta))/maxd))],1,d);
model.pd=pd; 
pr=randi([round((1-gamma)*((a2*landa*(1-betta))/maxr)), round((1+gamma)*((a2*landa*(1+betta))/maxr))],1,r);
model.pr=pr; 
ps=randi([round((1-gamma)*((a3*landa*(1-betta))/maxs)), round((1+gamma)*((a3*landa*(1+betta))/maxs))],1,s);
model.ps=ps;
pp=randi([round((1-gamma)*((sum(pn)*(1-betta))/maxp)), round((1+gamma)*((sum(pn)*(1+betta))/maxp))],1,p); 
model.pp=pp; 
pm=randi([round((1-gamma)*((sum(pp)*(1-betta))/m)), round((1+gamma)*((sum(pp)*(1+betta))/m))],1,m); 
model.pm=pm; 

% Opening cost of facilities 
betta=rand; 
gamma=rand; 
fcd=randi([round((1-betta)*(sum(sum(tcjd))+sum(sum(tcdp)))/maxd), round((1+betta)*(sum(sum(tcjd))+sum(sum(tcdp)))/maxd)], 1, d); 
model.fcd=fcd; 
fcm=randi([round((1-gamma)*(sum(sum(tcmp))+sum(sum(tcsm)))/maxm), round((1+gamma)*(sum(sum(tcmp))+sum(sum(tcsm)))/maxm)], 1, m); 
model.fcm=fcm; 
fcn=randi([round((1-betta)*(sum(sum(tcpn))+sum(sum(tcrn)))/maxn), round((1+betta)*(sum(sum(tcpn))+sum(sum(tcrn)))/maxn)], 1, n); 
model.fcn=fcn; 
fcp=randi([round((1-betta)*(sum(sum(tcpn))+sum(sum(tcmp)))/maxp), round((1+betta)*(sum(sum(tcpn))+sum(sum(tcmp)))/maxp)], 1, p); 
model.fcp=fcp; 
fcr=randi([round((1-gamma)*(sum(sum(tcrn))+sum(sum(tcjr)))/maxr), round((1+gamma)*(sum(sum(tcrn))+sum(sum(tcjr)))/maxr)], 1, r); 
model.fcr=fcr; 
fcs=randi([round((1-betta)*(sum(sum(tcsm))+sum(sum(tcjs)))/maxs), round((1+betta)*(sum(sum(tcsm))+sum(sum(tcjs)))/maxs)], 1, s);
model.fcs=fcs; 
fcj=randi([round((1-betta)*(sum(sum(tcjr))+sum(sum(tcjs)))/maxj), round((1+betta)*(sum(sum(tcjr))+sum(sum(tcjs)))/maxj)], 1, j);
model.fcj=fcj;

Qc=10000; % Financial risk 
model.Qc=Qc; 
end

