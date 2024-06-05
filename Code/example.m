%% Initialization
clc
clear 
warning off
close all

%% Data Processing
Nunits = 3;  %% Numbers of power generating plants
Horizon = 6; %% Time
Pmax = [100;50;25]; %% Maximum power capacity
Pmin = [20;10;0];   %% Minimum power capacity
C = [10 30 50];   %% Linear cost price
Pforecast = [80 130 160 150 120 80]; %% The forecasted power demand

%% OPF Model
x = sdpvar(Nunits,Horizon,'full');
Constraints = [];
for k = 1:Horizon
      Constraints = [Constraints, Pmin <= x(:,k) <= Pmax];
end
for k = 1:Horizon
      Constraints = [Constraints, sum(x(:,k)) == Pforecast(k)];
end
Objective = 0;
for k = 1:Horizon
  Objective = Objective + C*x(:,k);
end
ops = sdpsettings('solver', 'cplex');
optimize(Constraints, Objective, ops) % Set the solver to glpk

%% Visualization
bar(value(x)','stacked')
ylim([0,170])
ylabel('Power/MW')
xlabel('Time/h')
legend('Unit 1','Unit 2','Unit 3');
value(x)

%% Verification using MATPOWER (Optional)
% mpopt=mpoption;
% mpopt = mpoption(mpopt);
% rundcopf(mpc)