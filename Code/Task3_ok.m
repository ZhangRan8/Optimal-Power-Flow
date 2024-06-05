%% Initialization
clc
clear
warning off
close all

%% Data Processing
InputESS;
mpc = case39ee;
bus=mpc.bus;
gen=mpc.gen;
branch=mpc.branch;
gencost=mpc.gencost;
baseMVA=mpc.baseMVA;
UT=mpc.UT;
UD=mpc.UD;
Cs=mpc.Cs;

Nunits = 10;  %% Numbers of power generating plants
Nbus=39;
Horizon = 12; %% Time
Pmax = gen(:,9); %% Maximum power capacity
Pmin = gen(:,10);   %% Minimum power capacity
C = gencost(:,6)';   %% Linear cost price(c1)
Cl = gencost(:,7)';   %% c0
PLmax=branch(:,6);

% ----------------------------
% Pforecast = [80 130 160 150 120 80]; %% The forecasted power demand
LoadRate=ESS(:,1)';
Pforecast=bus(:,3)*LoadRate; % T=12 Load

WindRate=0.1;
TotalGen=WindRate*sum(Pmax);
Pw=ESS(:,2)*TotalGen;  % The predicted output of wind generator
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

% Form Ce (ESS matrix)
Ce=zeros(Nbus,1);
Ce(13,1)=1;


%% OPF Model
x = sdpvar(Nunits,Horizon,'full');
theta=sdpvar(Nbus,Horizon,'full');
xw= sdpvar(1,Horizon,'full');
% Unit binary variables
u=binvar(Nunits,Horizon,'full');
y=binvar(Nunits,Horizon,'full');
% ESS binary variables
xc=sdpvar(1,Horizon,'full'); % charge power of ESS at time
xs=sdpvar(1,Horizon,'full'); % discharge power of ESS at time
uc= binvar(1,Horizon,'full'); % the charging status of ESS
us= binvar(1,Horizon,'full'); %  the discharging status of ESS
E= sdpvar(1,Horizon,'full'); % the stored energy of ESS at time 

% Generator constraints
Constraints = [];
for k = 1:Horizon
    Constraints = [Constraints, Pmin .* u(:,k) <= x(:,k) <= Pmax.* u(:,k) ];
end

% Balancing constraints
balanceL=Bbus*theta;
balanceR=Cg*x + Cw*xw + Ce*(xs-xc)-Pforecast*LoadRatio;
Constraints = [Constraints, balanceL == balanceR];

% Line constraints
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

% Generator minimum up/down time constraints
for k = 2: Horizon
    Constraints = [Constraints, y(:,k) >= (u(:,k)-  u(:,k-1)) ];
end

for i = 1: Nunits
    t0=UT(i);
    for k = 2 : Horizon-t0 +1
        Constraints = [Constraints, sum(u(i,  k : k+t0-1)) >= t0*y(i,k)];
    end
    t0=UD(i);
    for k = 2 : Horizon-t0 +1
        Constraints = [Constraints, sum(u(i,  k : k+t0-1)) >= t0* ( u(i,k-1)-u(i,k) ) ];
    end
end

 Constraints = [Constraints, y(:,1)==0];  

%% ESS constraints
eps=0.1;
etac=0.9;
etas=0.9;
% Ecap=200;  
Ecap=200; 
E0=0.2*Ecap;


Psmax = 50;
Pcmax=Psmax;

for k=1:Horizon
    Constraints = [Constraints, uc(:, k)+us(:, k) <= 1];
    Constraints = [Constraints, 0<= xc(:,k) <= Pcmax * uc(:,k) ];
    Constraints = [Constraints, 0<= xs(:,k) <= Psmax * us(:,k) ];
    Constraints = [Constraints, eps*Ecap<= E(:,k) <= Ecap  ];  
end


for k=2:Horizon
    Constraints = [Constraints, E(:,k)==E(:,k-1) + etac*xc(:,k)-etas*xs(:,k) ];
end
Constraints = [Constraints, E(:,1)==E0 + etac*xc(:,1)-etas*xs(:,1) ];

%%
Objective = 0;
% Obj  f1
for k = 1:Horizon
    Objective = Objective + C*x(:,k) + Cl*u(:,k);
end
for k = 2:Horizon
    Objective = Objective + Cs*y(:,k);
end
% Objective = Objective +sum(Cl) * Horizon;

% Obj f2
Lambda=10;
f2=sum(Pw-xw);
f2=f2*Lambda;
Objective=Objective+f2;

ops = sdpsettings('solver', 'cplex');
optimize(Constraints, Objective, ops) % Set the solver to cplex
Totalvalue=value(Objective)

%% Calculate
% ObMy=value(Objective)
rl=Bf*value(theta)./(PLmax*LoadRatio);
rl=rl';

% result.x=value(x);
% result.u=value(u);
% result.y=value(y);
% 
% result.xc=value(xc);
% result.xs=value(xs);
% 
% result.uc=value(uc);
% result.us=value(us);
% 
% result.E=value(E);
%% Visualization

Ptotal=[value(x)];
Ptotal(Nunits+1,:)=value(xw);
% 创建堆叠图
b = bar(value(Ptotal)', 'stacked');
colors = parula(size(Ptotal, 2));
for i = 1:11
    set(b(i),'facecolor',colors(i,:))   
end


hold on
time=1:Horizon;
yP=sum(Pforecast)';
plot(time,yP,'LineWidth', 2)
legend('Unit 1','Unit 2','Unit 3','Unit 4','Unit 5','Unit 6','Unit 7','Unit 8','Unit 9','Unit 10','Wind','Load');

xlim([0,13])
ylim([0,8500])
ylabel('Power/MW')
xlabel('Time/h')


%% Calculate
SOC=value(E)/Ecap;
SOC=[0.2 SOC];
figure
time=1:Horizon+1;
plot(time,SOC)

%% Verification using MATPOWER (Optional)
% mpopt=mpoption;
% mpopt = mpoption(mpopt);
% rundcopf(mpc)