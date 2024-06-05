%% Initialization
clc
clear 
warning off
close all

%% Data Processing
InputWind;
mpc = case39ee;
bus=mpc.bus;
gen=mpc.gen;
branch=mpc.branch;
gencost=mpc.gencost;
baseMVA=mpc.baseMVA;

Nunits = 10;  %% Numbers of power generating plants
Nbus=39;
Horizon = 24; %% Time
Pmax = gen(:,9); %% Maximum power capacity
Pmin = gen(:,10);   %% Minimum power capacity
C = gencost(:,6)';   %% Linear cost price(c1)
Cl = gencost(:,7);   %% c0
PLmax=branch(:,6);

% ----------------------------
% Pforecast = [80 130 160 150 120 80]; %% The forecasted power demand
LoadRate=Wind(:,1)';
Pforecast=bus(:,3)*LoadRate; % T=24 Load

WindRate=0.5;
TotalGen=WindRate*sum(Pmax);
Pw=Wind(:,2)*TotalGen;  % The predicted output of wind generator
Pw=Pw';
% ------------------------------

% Set load ratio
LoadRatio=1;

% Form Bbus, Bf
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);

% Form Cg (Generator matrix)
col=gen(:,1);
Cg=zeros(Nbus,Nunits);
for i=1:Nunits
    Cg(col(i),i)=1;
end

% Form Cw (Wind generator matrix)
Cw=zeros(Nbus,1);
Cw(14,1)=1;

%% OPF Model
x = sdpvar(Nunits,Horizon,'full');
theta=sdpvar(Nbus,Horizon,'full');
xw= sdpvar(1,Horizon,'full');

% Generator constraints
Constraints = [];
for k = 1:Horizon
      Constraints = [Constraints, Pmin <= x(:,k) <= Pmax];
end

%Balancing constraints
balanceL=Bbus*theta;
balanceR=Cg*x + Cw*xw-Pforecast*LoadRatio;
Constraints = [Constraints, balanceL == balanceR];
for k = 1:Horizon
    temp = Bf*theta;
    Constraints = [Constraints, -PLmax*LoadRatio<=temp(:,k)<=PLmax*LoadRatio];
end

% slack bus
slack = find(bus(:,2)==3);
Constraints = [Constraints, theta(slack)==0];

% Ramping constraints
Eta=0.7;
RD=Eta*Pmax;
RU=RD;
for k = 2: Horizon
     Constraints = [Constraints, -RD<= (x(:,k)-  x(:,k-1)) <= RU ];
end

% Wind generator constraints
 Constraints = [Constraints, 0<= xw <= Pw ];

%%
Objective = 0;
% Obj  f1
for k = 1:Horizon
  Objective = Objective + C*x(:,k);
end
Objective = Objective +sum(Cl) * Horizon;
% Obj f2
Lambda=10;
f2=sum(Pw-xw);
f2=f2*Lambda;
Objective=Objective+f2;

ops = sdpsettings('solver', 'cplex');
optimize(Constraints, Objective, ops) % Set the solver to cplex

%% Calculate
% ObMy=value(Objective)
rl=Bf*value(theta)./(PLmax*LoadRatio);
rl=rl';

%% Visualization

Ptotal=[value(x)];
Ptotal(11,:)=value(xw);
% 创建堆叠图
b = bar(value(Ptotal)', 'stacked');
colors = parula(size(Ptotal, 2));
for i = 1:11
    set(b(i),'facecolor',colors(2*i,:))   
end


hold on
time=1:Horizon;
yP=sum(Pforecast)';
plot(time,yP,'LineWidth', 2)
legend('Unit 1','Unit 2','Unit 3','Unit 4','Unit 5','Unit 6','Unit 7','Unit 8','Unit 9','Unit 10','Wind','Load');

xlim([0,25])
ylim([0,8500])
ylabel('Power/MW')
xlabel('Time/h')

% legend('Load')
LoadRatio
x_re=value(x)



%% Verification using MATPOWER (Optional)
% mpopt=mpoption;
% mpopt = mpoption(mpopt);
% rundcopf(mpc)

