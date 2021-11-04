function dqdt=ode(t,QQ,pars,data,tov,v1tgb,v2tgb,v1apobb,v2apobb,pltg_b,basetg,mealtg,ss)

% function ode defines the system of dif. eq.

% pars: model parameters
% inl: injected labeled leucine per plasma vol.
% ing: injected labeled glycerol per plasma vol.
% data: patient data
% tov: time of visit, 1 for pre surgery and 2 for post surgery
% v1tgb, v2tgb, v1apobb, v2apobb steady state baseline values
% pltg: plasma tg
% pltg_b: steady state baseline plasma tg
% ss: the steady state constraint, 0 if no, 1 if yes
% basetg: baseline plasma tg
% mealtg: total amount of TG contained in the mixed meal


ins=QQ(80);
ins_b=data{2,12}(2,1);

cmtg_data=data{2,14}([1,2],:);
cmtg_b=cmtg_data(2,1);

cmtg=QQ(67);
vtg=sum(QQ([14,15,16,17]));
pltg=cmtg+vtg+basetg;

k=k_matrix(pars,QQ,ins,ins_b,v1tgb,v2tgb,v1apobb,v2apobb,pltg,pltg_b,ss);
% k is the matrix that holds transfer coefficients, where Q'=k*Q;

Q=QQ(1:28); % tracee
q=QQ(29:56); % tracer

dQ=k*Q;
dq=k*q;

Qdig=QQ(57:68); % digestive track
Qin=QQ(69:80); % insulin system

x=pars;
pdigest=[x(31),x(32)]; % digestive track parameters
kdig=kdigest(pdigest,tov); % digestive track coeficient matrix

plGly=x(21);
plLeu=x(22);

u1  = plLeu*(k(28,1)+k(2,1)*k(21,2)/(k(21,2)+k(1,2))); %constant leucine input that keeps Q1 at plLeu
u8 = plGly*(k(28,8)+k(10,8)); %constant glycerol input that keeps Q8 at plGly
dQ(1)=dQ(1)+u1;
dQ(8)=dQ(8)+u8;

%% insulin delay
% hmusig matrix holds the insulin input signal parameters estimated from
% the plasma insulin data, where time dependent insulin input signal is the
% bell shaped Gausian function (f) given below

hmusig=data{2,15};

alpha=hmusig(1);
beta=hmusig(2);
delta=hmusig(3);
ib=data{2,12}(2,1); % baseline insulin concentration

k1ins=2; % insulin transfer rate from the plasma compartment to the first delay compartment is fixed to the population average

f=@(x) alpha*exp(-((x-beta)/delta)^2); % insulin input signal estimated from plasma insulin data Eq. A10

dins(1)=f(t)-k1ins*(Qin(1)-ib); %Eq. A11

k2ins=x(29);
kv1=-k2ins*ones(11,1);
kv2=k2ins*ones(10,1);
km=diag(kv1)+diag(kv2,-1);
dins(2:12,1)=km*Qin(2:12);

dins(2)=dins(2)+k1ins*Qin(1);

%% Digestive track introducing post prandial switch
if tov==1;
    i3=0;
    i1=1;
    t_cm=15;% time it took for patient to to consume the meal in min before the surgery
else
    i3=1;
    i1=0;
    t_cm=30;% time it took for patient to to consume the meal in min after the surgery
end

k4=x(32); %CM_TG clearence rate

t_c=t_cm/60; % converting minutes to hours
tg_inr=mealtg/t_c; % average lipid input rate during the meal
ddig=kdig*Qdig;
ddig(1)=ddig(1)+i1*tg_inr;
ddig(3)=ddig(3)+i3*tg_inr;
ddig(11)=ddig(11)+k4*cmtg_b-k4*Qdig(11); % the term k4*cmtg_b ensures basal value to be equal to cmtg_b
%%
% dQ is the tracee system
% dq is the tracer system
% ddig is the digestive track system
% dins is the delyaed insulin system

dqdt=[dQ;dq;ddig;dins];
end

