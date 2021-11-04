function [sol1,sol2,sol3]=simul(pars,tspan,inl,ing,data,tov)

%% Setting initial conditions to steady state.
%All the regulatory pathways and pool sizes are constant during the
%preprandial state. Therefore, the initial states are calculated with
%linsolve. However, the linsolve function must be replaced by a nonlinear 
%solver if a nonlinear system is used for the preprandial state. This will
%increase the compuational cost dramaically.

% pars: model parameters
% ind: indices of estimated parameters (used if some parameters are fixed).
% inl: injected labeled leucine per plasma vol.
% ing: injected labeled glycerol per plasma vol.
% data: patient data
% tov: time of visit, 1 for pre surgery and 2 for post surgery

x=pars;

pw=(500/ing)*(885/92); 
% plasma volume is calculated from injected glycerol, where ing is the
% injected glycerol per liter plasma. a total of 500 mg labelled glycerol
% is injected for each patient.

plGly=x(21); % Steady state unlabelled plasma glycerol conc.
plLeu=x(22); % Steady state unlabelled plasma leucine conc.

ss=1; % ss=1 for the steady state calculations
k=k_matrix(pars,zeros(68,1),1,1,1,1,1,1,1,1,ss); 
% for the steady state calculations values of the state variables and
% their baseline values are set to 1, where the choice of 1 is arbitrary.
% Because, as long as the variable is equal to its baseline value it will 
% not impose a regulatory effect. During the preprandial state, the system
% is assumed to be in the steady-state

ac1=zeros(28,1);

u1  = plLeu*(k(28,1)+k(2,1)*k(21,2)/(k(21,2)+k(1,2))); %constant leucine input that keeps Q1 at plLeu
u8 = plGly*(k(28,8)+k(10,8)); %constant glycerol input that keeps Q8 at plGly

ac1(1)=u1;
ac1(8)=u8;
Qs1=linsolve(k,-ac1);

Q_ss=[Qs1;zeros(28,1);zeros(12,1);ones(12,1)]; % vector of initial states

pltg_b_d=data{2,11}(2,1); % baseline plasma tg

cmtg_data=data{2,14}([1,2],:);
cmtg_b_d=cmtg_data(2,1); % baseline plasma chylomicron tg

vtg_d=data{2,1}(2,:)+data{2,3}(2,:); % baseline total plasma vldl tg
vtg_b_d=mean(vtg_d(1:7));
basetg=pltg_b_d-vtg_b_d-cmtg_b_d; 
% basetg is the baseline plasma tg in the fractions other than 
% vldl and chylomicron. This is used to estimate total plasma 
% TG with the model.

v1tgb=sum(Q_ss([14,15])); % vldl1 tg baseline
v2tgb=sum(Q_ss([16,17])); % vldl2 tg baseline
v1apobb=sum(Q_ss([4,5])); % vldl1 apob baseline
v2apobb=sum(Q_ss([6,7])); % vldl2 apob baseline
vtg=sum(Q_ss([14,15,16,17])); % vldl tg baseline

pltg_b=cmtg_b_d+vtg+basetg; % plasma tg baseline

Q0=Q_ss(1:28); %initial values for the tracee system state variables
q0=zeros(28,1); % tracer values are 0 at t=0

Q0d=zeros(12,1);  % initial values for gastrointestinal module
Q0d(11)=cmtg_b_d; % intial value for plasma chylomicron Tg is set to the 
% baseline value obtained from the data

Q0in(1)=data{2,12}(2,1); % basal plasma insulin from data

k2ins=x(29);
Q0in(2:12,1)=(2/k2ins)*data{2,12}(2,1)*ones(11,1); % initial values for the 
%insulin module are set to steady state values, which is set equal to the basal
%insulin concentration for ease of use.

 
q0(1)=inl; % injected labelled leucine per liter plasma
q0(8)=ing; % injected labelled glycerol per liter plasma
Qq_0=[Q0;q0;Q0d;Q0in];

%%
ss=0; % ss is set to 0 for marking the nonsteady state calculations

if tov==1; % tov=1 for presurgery visit and tov=2 for post-surgery visit
    t_cm=15;% time that it took for patient to to consume the meal in min.
else
    t_cm=30;
end

t_c=t_cm/60;
mealtg=67.4*(100/95)*1000/pw; %TG from consumed lipids per plasma vol.;
% lipid content of the consumed meal: 67.4 gr
% 100 gr TG yields approximately 95 grams fatty acid (Ratnayake and Galli, Ann Nutr Metab 2009;55:8–43)


% ODE SOLUTION
% The system ODEs is solved in 3 steps;
% sol1: preprandial state
% sol2: meal consumption period
% sol3: postprandial state
% the end point of each set of solution is passed as initial condition to
% the next step to ensure continuity
opt=odeset('MaxStep',1);

sol1=ode15s(@(t,y) ode(t,y,pars,data,tov,v1tgb,v2tgb,v1apobb,v2apobb,pltg_b,basetg,0,ss), [tspan(1),2], Qq_0,opt);
Qq_0=sol1.y(:,end);
sol2=ode15s(@(t,y) ode(t,y,pars,data,tov,v1tgb,v2tgb,v1apobb,v2apobb,pltg_b,basetg,mealtg,ss), [2,2+t_c], Qq_0,opt);
Qq_0=sol2.y(:,end);
sol3=ode15s(@(t,y) ode(t,y,pars,data,tov,v1tgb,v2tgb,v1apobb,v2apobb,pltg_b,basetg,0,ss), [2+t_c,tspan(2)], Qq_0,opt);

end