function fig_pl(par_est,tspan,inl,ing,datameans,tov)

x=par_est;

%% simulating the model with given parameters
[sol1,sol2,sol3]=simul(par_est,tspan,inl,ing,datameans,tov);

%% Getting Data
t1=datameans{2,1}(1,:);

v1tgp=datameans{2,1}(2,:);
v2tgp=datameans{2,3}(2,:);
v1tgen=datameans{2,2}(2,:);
v2tgen=datameans{2,4}(2,:);

v1apop=datameans{2,5}(2,:);
v2apop=datameans{2,7}(2,:);
v1apopen=datameans{2,6}(2,:);
v2apopen=datameans{2,8}(2,:);

plleuen=datameans{2,9}(2,:);
plglyen=datameans{2,10}(2,:);

v1tgpe=datameans{2,1}(3,:);
v2tgpe=datameans{2,3}(3,:);
v1tgene=datameans{2,2}(3,:);
v2tgene=datameans{2,4}(3,:);

v1apope=datameans{2,5}(3,:);
v2apope=datameans{2,7}(3,:);
v1apopene=datameans{2,6}(3,:);
v2apopene=datameans{2,8}(3,:);

plleuene=datameans{2,9}(3,:);
plglyene=datameans{2,10}(3,:);

tpltg=datameans{2,11}(1,:);
pltg=datameans{2,11}(2,:);
spltg=datameans{2,11}(3,:);

cmtg=datameans{2,14}(2,:);
scmtg=datameans{2,14}(2,:);

vtg=v1tgp+v2tgp;
vtg0=mean(vtg(1:7));
basetg=pltg(1)-vtg0-cmtg(1);

%basetg=x(33);

%cmtg=pltg-interp1(vtg,t1,tpltg)-basetg; %actual CM_TG data is used

%% evaluating solution at data points;
%There are 4 different measurement time sets;
% t1 is fot vltgp, v2tgp, v1apop, v2apop
% t2 is for v1tgen and v2tgen
% t3 is for plleuen
% t4 is for plglyen
% t5 is for v1apopen, v2apopen
% t6 is for pltg


t1=datameans{2,1}(1,:);
t2=datameans{2,2}(1,:);
t3=datameans{2,9}(1,:);
t4=datameans{2,10}(1,:);
t5=datameans{2,6}(1,:);
t6=datameans{2,11}(1,:);
t7=datameans{2,14}(1,:);

tsim=tspan(1):0.01:tspan(2);
Q=sols(sol1,sol2,sol3,tsim);

ci=[4,5];
v1p=ci;
v2p=ci+2;
v1t=ci+10;
v2t=ci+12;

v1apop_s    =sum(Q(v1p,:));
v2apop_s    =sum(Q(v2p,:));
v1tgp_s     =sum(Q(v1t,:));
v2tgp_s     =sum(Q(v2t,:));

vtgp_s=sum(Q(v1t,:))+sum(Q(v2t,:));
cmtg_s=Q(67,:);
pltg_s=vtgp_s+cmtg_s+basetg;

nc=28;

v1apopen_s  =sum(Q(v1p+nc,:))./(sum(Q(v1p+nc,:))+sum(Q(v1p,:)))*100; %TTR
v2apopen_s  =sum(Q(v2p+nc,:))./(sum(Q(v2p+nc,:))+sum(Q(v2p,:)))*100; %TTR

v1tgen_s    =sum(Q(v1t+nc,:))./(sum(Q(v1t+nc,:))+sum(Q(v1t,:)))*100; %APE data is not multiplied by 100
v2tgen_s    =sum(Q(v2t+nc,:))./(sum(Q(v2t+nc,:))+sum(Q(v2t,:)))*100; %APE data is not multiplied by 100

plleuen_s   =100*Q(1+nc,:)./(Q(1+nc,:)+Q(1,:));     %APE
plglyen_s   =Q(8+nc,:)./(Q(8,:));                    %TTR

vtgp=v1tgp+v2tgp;
vapop=v1apop+v2apop;
vtgen=v1tgen+v2tgen;
vapopen=v1apopen+v2apopen;

vtgpe=v1tgpe;
vapope=v1apope;
vtgene=v1tgene;
vapopene=v1apopene;

vtgp_s=v1tgp_s+v2tgp_s;
vapop_s=v1apop_s+v2apop_s;
vtgen_s=v1tgen_s+v2tgen_s;
vapopen_s=v1apopen_s+v2apopen_s;

%%
figure(2)
subplot(3,2,1)
plot(tsim,v1apopen_s)
hold on;
errorbar(t5,v1apopen,v1apopene,'*')
title('VLDL1 ApoB APE')

subplot(3,2,2)
plot(tsim,v2apopen_s)
hold on;
errorbar(t5,v2apopen,v2apopene,'*')
title('VLDL2 ApoB APE')

subplot(3,2,3)
plot(tsim,v1tgen_s)
hold on;
errorbar(t2,v1tgen,v1tgene,'*')
title('VLDL1 TG APE')

subplot(3,2,4)
plot(tsim,v2tgen_s)
hold on;
errorbar(t2,v2tgen,v2tgene,'*')
title('VLDL2 TG APE')

subplot(3,2,5)
plot(tsim,plleuen_s)
hold on;
errorbar(t3,plleuen,plleuene,'*')
title('PL LEU APE')

subplot(3,2,6)
plot(tsim,plglyen_s)
hold on;
errorbar(t4,plglyen,plglyene,'*')
title('PL GLY APE')


figure(1)

subplot(3,2,1)
plot(tsim,v1apop_s)
hold on;
errorbar(t1,v1apop,v1apope,'*')
title('VLDL1 ApoB mg/l')
ylabel('VLDL_1 APoB mg/L');

subplot(3,2,2)
plot(tsim,v2apop_s)
hold on;
errorbar(t1,v2apop,v2apope,'*')
title('VLDL2 ApoB mg/l')
ylabel('VLDL_2 APoB mg/L');

subplot(3,2,3)
plot(tsim,v1tgp_s)
hold on;
errorbar(t1,v1tgp,v1tgpe,'*')
title('VLDL1 TG mg/l')
ylabel('VLDL_1 TG mg/L');

subplot(3,2,4)
plot(tsim,v2tgp_s)
hold on;
errorbar(t1,v2tgp,v2tgpe,'*')
title('VLDL2 TG mg/l')
ylabel('VLDL_2 TG mg/L');

subplot(3,2,5)
plot(tsim,cmtg_s)
hold on;
errorbar(t7,cmtg,scmtg,'*')
title('CM TG')
ylabel('CM TG mg/L');

subplot(3,2,6)
plot(tsim,pltg_s)
hold on;
errorbar(t6,pltg,spltg,'*')
title('Total Plasma TG')
ylabel('Total Plasma TG mg/L');


end

