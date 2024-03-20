function z= MyCost(x,model)
%% Input indices and parameters 
m=model.m; 
p=model.p; 
n=model.n; 
i=model.i; 
d=model.d; 
r=model.r; 
j=model.j; 
s=model.s; 
c=model.c; 
pbc=model.pbc;
Qc=model.Qc; 
accni=model.accni;
mcpc=model.mcpc;
pcmc=model.pcmc; 
tcmp=model.tcmp; 
tcpn=model.tcpn;
tcsm=model.tcsm;
tcjs=model.tcjs; 
tcjd=model.tcjd;
tcdp=model.tcdp;
tcjr=model.tcjr;
tcrn=model.tcrn;
tcij=model.tcij;
maxd=model.maxd;
maxm=model.maxm;
maxn=model.maxn;
maxp=model.maxp;
maxr=model.maxr; 
maxs=model.maxs;
maxj=model.maxj; 
dci=model.dci;
bp=model.bp;
alphai=model.alphai;
scj=model.scj;
scd=model.scd; 
scr=model.scr; 
a1=model.a1; 
a2=model.a2;
a3=model.a3; 
pn=model.pn; 
pj=model.pj;
pd=model.pd;
pr=model.pr;
ps=model.ps;
ppp=model.pp; 
pm=model.pm; 
fcd=model.fcd;
fcm=model.fcm;
fcn=model.fcn; 
fcp=model.fcp; 
fcr=model.fcr; 
fcs=model.fcs;
fcj=model.fcj;
%% Variables 

%Location decisions 
Yd=zeros(1,d); 
Xd=x(1:d); 
[a, b]=sort(Xd); 
for dd=1:d
    if b(dd)<=maxd
        Yd(dd)=1;
    else
        Yd(dd)=0; 
    end
end
a=[]; 
b=[]; 
Ym=zeros(1, m); 
Xm=x(d+1: d+m); 
[a, b]=sort(Xm); 
for mm=1:m
    if b(mm)<=maxm
        Ym(mm)=1;
    else
        Ym(mm)=0; 
    end
end
a=[]; 
b=[]; 
Yn=zeros(1, n); 
Xn=x(d+m+1: d+m+n); 
[a, b]=sort(Xn);
for nn=1:n
    if b(nn)<=maxn
        Yn(nn)=1; 
    else
        Yn(nn)=0; 
    end
end
Yp=zeros(1,p); 
Xp=x(d+m+n+1: d+m+n+p); 
[a, b]=sort(Xp); 
for pp=1:p
    if b(pp)<= maxp
        Yp(pp)=1; 
    else
        Yp(pp)=0; 
    end
end
Yr=zeros(1,r); 
Xr=x(d+m+n+p+1: d+m+n+p+r); 
[a, b]=sort(Xr);
for rr=1:r
    if b(rr)<=maxr
        Yr(rr)=1;
    else
        Yr(rr)=0; 
    end
end
Ys=zeros(1,s);
Xs=x(d+m+n+p+r+1: d+m+n+p+r+s); 
[a, b]=sort(Xs); 
for ss=1:s
    if b(ss)<=maxs
        Ys(ss)=1; 
    else
        Ys(ss)=0; 
    end
end
Yj=zeros(1, j); 
Xj=x(d+m+n+p+r++s+1: d+m+n+p+r+s+j); 
[a, b]=sort(Xj); 
for jj=1:j
    if b(jj)<=maxj
        Yj(jj)=1; 
    else
        Yj(jj)=0; 
    end
end

% Allocation decisions 
Xijc=zeros(i,j,c); 
Xjrc=zeros(j,r,c); 
Xjdc=zeros(j,d,c); 
Xjsc=zeros(j,s,c); 
Xsmc=zeros(s,m,c); 
Xmpc=zeros(m,p,c); 
Xdpc=zeros(d,p,c); 
Xpnc=zeros(p,n,c); 
Xrnc=zeros(r,n,c); 
Xnic=zeros(n,i,c); 
for ii=1:i
    for nn=1:n
        for cc=1:c
            if Yn(nn)==1
                if pn(nn)> dci(cc, ii)
                    Xnic(nn,ii,cc)=Xnic(nn,ii,cc)+dci(cc,ii); 
                    pn(nn)=pn(nn)-dci(cc,ii); 
                    dci(cc,ii)=0; 
                end
            end
        end
    end
end
for nn=1:n
    for pp=1:p
        for cc=1:c
            if Yp(pp)==1 && Yn(nn)==1
                if ppp(pp)>sum(Xnic(nn,:,cc))
                    Xpnc(pp,nn,cc)=sum(Xnic(nn,:,cc))+Xpnc(pp,nn,cc); 
                    ppp(pp)=ppp(pp)-sum(Xnic(nn,:,cc)); 
                end
            end
        end
    end
end
for mm=1:m
    for pp=1:p
        for cc=1:c
            if Ym(mm)==1 && Yp(pp)==1
                if pm(mm)>sum(Xpnc(pp,:,cc))
                    Xmpc(mm,pp,cc)=Xmpc(mm,pp,cc)+sum(Xpnc(pp,:,cc)); 
                    pm(mm)=pm(mm)-sum(Xpnc(pp,:,cc)); 
                end
            end
        end
    end
end
for ii=1:i
    for jj=1:j
        for cc=1:c
            if Yj==1 
                if pj(jj)> (model.dci(ii,cc))*alphai(ii)
                    Xijc(ii,jj,cc)=Xijc(ii,jj,cc)+(model.dci(ii,cc))*alphai(ii); 
                    pj(jj)=pj(jj)-(model.dci(ii,cc))*alphai(ii); 
                    alphai(ii)=0; 
                end
            end
        end
    end
end

for jj=1:j
    for rr=1:r
        for cc=1:c
            if Yj(jj)==1 && Yr(rr)==1
                if pr(rr)>sum(Xijc(:,jj,cc))*a1
                    Xjrc(jj,rr,cc)=Xjrc(jj,rr,cc)+sum(Xijc(:,jj,cc))*a1; 
                    pr(rr)=pr(rr)-sum(Xijc(:,jj,cc))*a1; 
                end
            end
        end
    end
end
for jj=1:j
    for dd=1:d
        for cc=1:c
            if Yj(jj)==1 && Yd(dd)==1
                if pd(dd) >sum(Xijc(:,jj,cc))*a2
                    Xjdc(jj,dd,cc)=Xjdc(jj,dd,cc)+sum(Xijc(:,jj,cc))*a2; 
                    pd(dd)=pd(dd)-sum(Xijc(:,jj,cc))*a2; 
                end
            end
        end
    end
end
for jj=1:j
    for ss=1:s
        for cc=1:c
            if Yj(jj)==1 && Ys(ss)==1
                if ps(ss)> sum(Xijc(:,jj,cc))*a3
                    Xjsc(jj,ss,cc)=Xjsc(jj,ss,cc)+sum(Xijc(:,jj,cc))*a3; 
                    ps(ss)=ps(ss)-sum(Xijc(:,jj,cc))*a3; 
                end
            end
        end
    end
end
for ss=1:s
    for mm=1:m
        for cc=1:c
            if Ys(ss)==1 && Ym(mm)==1
                Xsmc(ss,mm,cc)=sum(Xjsc(:,ss,cc)); 
            end
        end
    end
end
for dd=1:d
    for pp=1:p
        for cc=1:c
            if Yd(dd)==1 && Yp(pp)==1
                Xdpc(dd,pp,cc)=Xdpc(dd,pp,cc)+sum(Xjdc(:,dd,cc)); 
            end
        end
    end
end
for rr=1:r
    for nn=1:n
        for cc=1:c
            if Yr(rr)==1 && Yn(nn)==1
                Xrnc(rr,nn,cc)=Xrnc(rr,nn,cc)+sum(Xjrc(:,rr,cc)); 
            end
        end
    end
end

%% Objective functions 
z1=0; 
for dd=1:d
    z1=z1+Yd(dd)*fcd(dd); 
end
for mm=1:m
    z1=z1+Ym(mm)*fcm(mm); 
end
for nn=1:n
    z1=z1+Yn(nn)*fcn(nn); 
end
for pp=1:p
    z1=z1+fcp(pp)*Yp(pp); 
end
for rr=1:r
    z1=z1+fcr(rr)*Yr(rr); 
end
for ss=1:s
    z1=z1+fcs(ss)*Ys(ss); 
end
for jj=1:j
    z1=z1+fcj(jj)*Yj(jj); 
end
TCAC=zeros(1,c); 
for cc=1:c
    for nn=1:n
        for ii=1:i
            TCAC(cc)=TCAC(cc)+Xnic(nn,ii,cc)*accni(nn,ii,cc); 
        end
    end
end
TCMC=zeros(1,c);
for cc=1:c
    for pp=1:p
        TCMC(cc)=TCMC(cc)+sum(Xpnc(pp,:,cc))*mcpc(pp,cc)*bp(pp); 
    end
end
TCPC=zeros(1,c); 
for cc=1:c
    for pp=1:p
        TCPC(cc)=TCPC(cc)+sum(Xpnc(pp,:,cc))*pcmc(pp,cc); 
    end
end
TCTC=zeros(1,c); 
for cc=1:c
    for pp=1:p
        for nn=1:n
            TCTC(cc)=TCTC(cc)+Xpnc(pp,nn,cc)*tcpn(pp,nn); 
        end
    end
end
for cc=1:c
    for pp=1:p
        for mm=1:m
            TCTC(cc)=TCTC(cc)+Xmpc(mm,pp,cc)*tcmp(mm,pp); 
        end
    end
end
for cc=1:c
    for ii=1:i
        for jj=1:j
            TCTC(cc)=TCTC(cc)+tcij(ii, jj)*Xijc(ii,jj,cc); 
        end
    end
end
for cc=1:c
    for jj=1:j
        for rr=1:r
            TCTC(cc)=TCTC(cc)+Xjrc(jj,rr,cc)*tcjr(jj,rr); 
        end
    end
end
for cc=1:c
    for jj=1:j
        for dd=1:d
            TCTC(cc)=TCTC(cc)+Xjdc(jj,dd,cc)*tcjd(jj,dd); 
        end
    end
end
for cc=1:c
    for jj=1:j
        for ss=1:s
            TCTC(cc)=TCTC(cc)+Xjsc(jj,ss,cc)*tcjs(jj,ss);
        end
    end
end
for cc=1:c
    for ss=1:s
        for mm=1:m
            TCTC(cc)=TCTC(cc)+Xsmc(ss,mm,cc)*tcsm(ss,mm); 
        end
    end
end
for cc=1:c
    for dd=1:d
    for pp=1:p
        TCTC(cc)=TCTC(cc)+tcdp(dd,pp)*Xdpc(dd,pp,cc); 
    end
    end
end
for cc=1:c
    for rr=1:r
        for nn=1:n
           TCTC(cc)=TCTC(cc)+Xrnc(rr,nn,cc)*tcrn(rr,nn); 
        end
    end
end 
PRC=zeros(1,c); 
for cc=1:c
    for jj=1:j 
        PRC(cc)=PRC(cc)+scj(jj)*sum(Xijc(:,jj,cc)); 
    end
end
for cc=1:c
    for dd=1:d
        PRC(cc)=PRC(cc)+scd(dd)*sum(Xjdc(:,dd,cc)); 
    end
end
for cc=1:c
    for rr=1:r
        PRC(cc)=PRC(cc)+scr(rr)*sum(Xjrc(:,rr,cc)); 
    end
end

z11=zeros(1,c); 
for cc=1:c
    z11(cc)=pbc(cc)*(TCAC(cc)+TCMC(cc)+TCPC(cc)+TCTC(cc)-PRC(cc));
end
z1=sum(z11)+z1; 
z2=0; 
% for cc=1:c
%     if z11(cc)>Qc
%     z2=z2+(z11(cc)-Qc); 
%     end
% end
z=[z1, z2]';
end

