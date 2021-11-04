function k=kdigest(p,tov)
% Function creates the matrix that holds the transfer coeffcients for the
% gastrointestinal module.
% p is the parameter vector
% tov is the time of visit: 1 for pre-surgery and 2 for post-surgery


if tov==1;
    ma=0.08; % ratio of fat in the feces. (malabsorbtion coef. Odstrcil et al. 2010)
    dc=6;   % number of intestinal compartments
else
    ma=0.27; % ratio of fat in the feces. (malabsorbtion coef.) is raised to 0.3 after the surgery (Odstrcil et al. 2010)
    dc=4; % number of intestinal compartments are reduced to 4 after the surgery
end

rt=exp(log(ma)/dc);

k1=2.005;
k2=p(1);
k3=k2*(1-rt)/rt;
k4=p(2);

k=zeros(12,12);

for i=1:3; % Stomach compartments
    k(i+1,i)=k1;
end

for i=4:9; % intestinal compartments
    k(i+1,i)=k2;
end

% in the post-surgical situation last gastric compartment is directly 
% attached to the second intestinal compartment

if tov==2; 
    k(4,3)=0; 
    k(5,3)=k1;
end

k(11,[4,5,6,7,8,9])=k3; % intestinal lipid absorbtion

% intestinal lipid absorption from the first intestinal copartment after
% the surgery is set to 0 becasue when left as a free parameter its value
% approached to 0.

if tov==2;
    k(11,5)=0;
end

% compartment 12 is the exit compartment of 
% the gastrointestnal tract module. 

k(12,11)=k4; 
a=sum(k);

for i=1:12;
    k(i,i)=-a(i);
end

k(11,11)=0;
k(10,10)=0;
k(12,12)=-1; % this value is arbirarily set to -1 to make matrix non-singular and has no impact.
end