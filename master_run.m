% This file can be used to run the model
tspan=[0,10];
tov=[1,2]; % time of visit. 1 fror first visit (pre-surgery) and 2 for second visit (post-surgery)

load par_table.mat; 
% Tp is the table for parameter values. (Parameter values estimated for the population average)
% lb,ub; lower and upper bounds; 
% est1,est2 pre and post surgery estimated values; 
% std1 and std2 estimated parameter std.
% r1 and r1 are std/est ratios for pre- and post-surgery, respectively

load mdata.mat; % data for population averages;
load inp.mat; % injected labelled materials;

for i=tov;
    pars=Tp{:,5+(i-1)*2};
    inl=inp{1+i,1};
    ing=inp{1+i,2};
    datameans=mdata{1,i}; 
    fig_pl(pars,tspan,inl,ing,datameans,i);
end