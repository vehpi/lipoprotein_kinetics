function k=k_matrix(par_est,QQ,ins,ins_b,v1tgb,v2tgb,v1apobb,v2apobb,pltg,pltg_b,ss);

% The function creates the k matrix that holds the transfer coefficients 
% for the liver, tracer injection and plasma modules.
% par_est: model parameters
% QQ: sate variables (compartments 1:28 tracee, compartments 29:56 tracer)
% ins: plasma insulin 
% ins_b: baseline insulin
% v1tgb, v2tgb, v1apobb, v2apobb steady state baseline values
% pltg: plasma tg
% pltg_b: baseline total plasma tg
% ss: the steady state , 0 for non-steady state, 1 for steady state

num_comp=28; % The number of compartments in the hepatic, plasma and insulin injection modules
% The compartment 28 is the exit compartment

Q=QQ(1:28); % Tracee system compartments
k=zeros(num_comp);
x=par_est;

%% Competition for lipolysis between lipoproteins
n_inh=x(30);

v1tg=Q(14)+Q(15);
v2tg=Q(16)+Q(17);
v1apob=Q(4)+Q(5);
v2apob=Q(6)+Q(7);

rtg=pltg/pltg_b;

sv1b=(v1tgb)/v1apobb;
sv2b=(v2tgb)/v2apobb;

sv1=(v1tg)/v1apob;
sv2=(v2tg)/v2apob;

if ss==1;
    inh1=1;
    inh2=1;
else
    inh1=((sv1/rtg)*(1/sv1b))^n_inh;
    inh2=((sv2/rtg)*(1/sv2b))^n_inh;
end
%%

%parameters fitted during pre-prandial state;
d4_27=x(1);
delayTime=x(2);
delayTimeTG=x(3);
f5_4=x(4);
f6_4=x(5);
f7_6=x(6);
k(28,1)=x(7);
k(28,4)=x(8);
k(28,5)=x(9);
k(28,6)=x(10);
k(1,2)=x(11);
k(1,3)=x(12);
k(3,1)=x(13);
k(5,4)=x(14);
k(6,4)=x(15);
k(12,10)=x(16);
k(10,8)=x(17);
k(10,11)=x(18);
k(11,10)=x(19);
k(21,2)=x(20);

plasmaGlycerol=x(21);
plasmaLeucine=x(22);
kins=x(23);

k(28,8)=x(24);
k(9,8)=x(25);
k(8,9)=x(26);

k(7,6)=x(27);
k(28,7)=x(28);

k2ins=x(29);

% TG secretion per ApoB particle from compartment 20 to 16 is set to be
% equal to the TG per ApoB secretion from compartment 14 to 16.
d14_20=d4_27/(d4_27+(1-d4_27)*f6_4);

kD=7/delayTime; % leucine delay time constant
kDTG=5/delayTimeTG; % glycerol delay time constant

% leucine equations;

k(2,1)=k(1,2);

% Leucine delay
k(4,27)=d4_27*kD;
k(6,27)=(1-d4_27)*kD;

k(22,21)=kD;
k(23,22)=kD;
k(24,23)=kD;
k(25,24)=kD;
k(26,25)=kD;
k(27,26)=kD;

% Glycerol delay
k(14,20)=d14_20*kDTG;
k(16,20)=(1-d14_20)*kDTG;

k(13,12)=kDTG;
k(18,13)=kDTG;
k(19,18)=kDTG;
k(20,19)=kDTG;

%% Incorporating the competition in lipolysis and insulin mediated stimulation 

if ss==1;
    ins_b_ef=ins_b;
else
    ins_b_ef=(2/k2ins)*ins_b;
end

nins=1;

Rmax=0.25; % maximum LPL stimulation by insulin is set to 0.25
lpla=1+Rmax*((ins-ins_b_ef)^nins)/((ins-ins_b_ef)^nins+kins^nins); % LPL stimulation by insulin

k(28,6)=k(28,6)*inh2*lpla;
k(28,7)=k(28,7)*inh2*lpla;
k(7,6)=k(7,6)*inh2*lpla;

k(6,4)=k(6,4)*inh1*lpla;
k(28,4)=k(28,4)*inh1*lpla;
k(5,4)=k(5,4)*inh1*lpla;
k(28,5)=k(28,5)*inh1*lpla;

%%
% glycerol equations;
% glycerol (TG) system is coupled to the leucine (apoB) system via transfer
% rates
k(28,14)=(1-f6_4)*k(6,4) + k(28,4) + (1-f5_4)*k(5,4);
k(28,16)=k(28,6) + (1-f7_6)*k(7,6);

k(16,14)=k(6,4)*f6_4;

k(28,15)=k(28,5);
k(28,17)=k(28,7);
k(15,14)=k(5,4)*f5_4;
k(17,16)=k(7,6)*f7_6;

%% loss of material from compartments.
% setting the diagonal of the matrix to the sum of corresponding column 
% gives the loss of material from the corresponding compartment. 

a=sum(k);
k=k-diag(a);
k(28,28)=-1; % Kompartment 28 is the exit compartment, the arbitrary non-zero value is assigned to make matrix k non-singular and has no impact

end
