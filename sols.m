function y=sols(sol1,sol2,sol3,t);
% The function takes seperate consequtive sol structures obtained from 
% the simul function and concatenates them into a continuous data array y.

ts1=[sol1.x(1),sol1.x(end)];
ts2=[sol2.x(1),sol2.x(end)];
ts3=[sol3.x(1),sol3.x(end)];

y1=deval(sol1,t(ts1(1)<=t   & t<ts1(2)));
y2=deval(sol2,t(ts2(1)<=t   & t<=ts2(2)));
y3=deval(sol3,t(ts3(1)<t    & t<=ts3(2)));

y=[y1,y2,y3];

end
