%% Initialization
clc
clear 
warning off
close all

%% Data Processing
mpc = case39ee;
bus=mpc.bus;
gen=mpc.gen;
branch=mpc.branch;
gencost=mpc.gencost;
baseMVA=mpc.baseMVA;

Nunits = 10;  %% Numbers of power generating plants
Nbus=39;
Horizon = 1; %% Time (Single time)
Pmax = gen(:,9); %% Maximum power capacity
Pmin = gen(:,10);   %% Minimum power capacity
C = gencost(:,6)';   %% Linear cost price(c1)
Cl = gencost(:,7);   %% c0
PLmax=branch(:,6);
% Pforecast = [80 130 160 150 120 80]; %% The forecasted power demand
Pforecast=bus(:,3); %single time OPF only has one load
% Set load ratio (Task C)
LoadRatio=0.8;

% Form Bbus, Bf
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);

% Form Cg (Generator matrix)
col=gen(:,1);
Cg=zeros(Nbus,Nunits);
for i=1:Nunits
    Cg(col(i),i)=1;
end


%% OPF Model
x = sdpvar(Nunits,Horizon,'full');
theta=sdpvar(Nbus,Horizon,'full');
Constraints = [];
for k = 1:Horizon
      Constraints = [Constraints, Pmin <= x(:,k) <= Pmax];
end
% for k = 1:Horizon
%       Constraints = [Constraints, sum(x(:,k)) == Pforecast(k)];
% end


balanceL=Bbus*theta;
balanceR=Cg*x-Pforecast*LoadRatio;
Constraints = [Constraints, balanceL == balanceR];
Constraints = [Constraints, -PLmax*LoadRatio<=Bf*theta<=PLmax*LoadRatio];
slack = find(bus(:,2)==3);
Constraints = [Constraints, theta(slack)==0];


Objective = 0;
for k = 1:Horizon
  Objective = Objective + C*x(:,k);
end
Objective = Objective +sum(Cl);
ops = sdpsettings('solver', 'cplex');
optimize(Constraints, Objective, ops) % Set the solver to cplex

%% Calculate
% ObMy=value(Objective)
rl=Bf*value(theta)./(PLmax*LoadRatio)


%% Visualization
b=bar(value(x)','stacked');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0,1300])
ylabel('Power/MW')
xlabel('Unit')
% legend('Unit 1','Unit 2','Unit 3');
LoadRatio
%%
x_result=value(x);
theta_result=value(theta);
PL_result=Bf*theta_result;

%% Verification using MATPOWER (Optional)
mpopt=mpoption;
mpopt = mpoption(mpopt);
result=rundcopf(mpc)

%%
theta_deg = roundn(value(theta_result)*180/pi,-3)

