% MPC Code
% Jinsung Kim

clear; clc;

% Process parameters
A1=28; % cm
A2=32;
A3=28;
A4=32;

a1=0.071; % cm
a2=0.057;
a3=0.071;
a4=0.057;

kc=0.5; % V/cm
g=981;  % cm/s^2

v10_nmp=3.00; % V
v20_nmp=3.00;

k1_nmp=3.33; % cm^3/Vs
k2_nmp=3.35; % cm^3/Vs

g1_nmp=0.25; % ratio of allocated punmp capacity between lower and
            % upper tank
g2_nmp=0.35;

% If symbolic toolbox is installed, this code could be used to
% calculate the stationary values
% Stationary values, non-minimum phase
%S = solve('-a1/A1*sqrt(2*g*h10_nmp)+a3/A1*sqrt(2*g*h30_nmp)+g1_nmp*k1_nmp/A1*v10_nmp',...
%	  '-a2/A2*sqrt(2*g*h20_nmp)+a4/A2*sqrt(2*g*h40_nmp)+g2_nmp*k2_nmp/A2*v20_nmp',...
%	  '-a3/A3*sqrt(2*g*h30_nmp)+(1-g2_nmp)*k2_nmp/A3*v20_nmp',...
%	  '-a4/A4*sqrt(2*g*h40_nmp)+(1-g1_nmp)*k1_nmp/A4*v10_nmp',...
%	  'h10_nmp,h20_nmp,h30_nmp,h40_nmp');

%h10_nmp = eval(S.h10_nmp);
%h20_nmp = eval(S.h20_nmp);
%h30_nmp = eval(S.h30_nmp);
%h40_nmp = eval(S.h40_nmp);

h10_nmp = 8.24441415257276;
h20_nmp = 19.01629576919927;
h30_nmp = 4.31462580236556;
h40_nmp = 8.80652939083585;


% Build state space model, minimum phase
T1_nmp=A1/a1*sqrt(2*h10_nmp/g);
T2_nmp=A2/a2*sqrt(2*h20_nmp/g);
T3_nmp=A3/a3*sqrt(2*h30_nmp/g);
T4_nmp=A4/a4*sqrt(2*h40_nmp/g);

A_nmp=[-1/T1_nmp 0 A3/(A1*T3_nmp) 0;
      0 -1/T2_nmp 0 A4/(A2*T4_nmp);
      0 0 -1/T3_nmp 0;
      0 0 0 -1/T4_nmp];

B_nmp=[g1_nmp*k1_nmp/A1 0;
      0 g2_nmp*k2_nmp/A2;
      0 (1-g2_nmp)*k2_nmp/A3;
      (1-g1_nmp)*k1_nmp/A4 0];

C_nmp=[kc 0 0 0; % Notice the measured signals are given in Volts!
      0 kc 0 0];

D_nmp=zeros(2,2);

Ts = 4;

% Constraints
% No constraints on du
% Pump capacities [0 10]V
% Level 1 [0 20]cm = [0 10]V
% Level 2 [0 20]cm
% Level 3 [0 20]cm
% Level 4 [0 20]cm

% du_max = [10000 10000]'; % limit on delta u; slew rate
% du_min = [-10000 -10000]';
du_max = [1000 1000]'; % limit on delta u; slew rate
du_min = [-1000 -1000]';
u_max = [10-v10_nmp 10-v20_nmp]'; % limit absolute value of u
u_min = [-v10_nmp -v20_nmp]';
z_max = [10-h10_nmp/2-0.1 10-h20_nmp/2-0.1 10-h30_nmp/2-0.1 10-h40_nmp/2-0.1]'; % Limits on controlled outputs
z_min = [-h10_nmp/2 -h20_nmp/2 -h30_nmp/2 -h40_nmp/2]'; 

% Set point trajectory
ref = [zeros(60/Ts,1); 
    3*ones(1400/Ts,1)];
ref = [ref zeros(length(ref),1)];

% Input disturbance trajectory
disturbance = [zeros(600/Ts,1); -1*ones(860/Ts,1)];
disturbance = [zeros(length(disturbance),1) disturbance];

% ----  Control Setting ----
% MPC parameters
Np = 30; % Prediction horizon
Nu = 10; % Horizon for varying input signal

Q = diag([4 1]);
R = 0.01*diag([1 1]);

QQ = kron(eye(Np), Q);
RR = kron(eye(Nu), R);

[Ad, Bd, Cyd, Dzd]=ssdata(c2d(ss(A_nmp,B_nmp,C_nmp,D_nmp), Ts));
Bwd = zeros(size(Bd));
[nx, nu] = size(Bd);
Czd = Cyd;
nz = size(Czd, 1);

Ccd = 0.5*eye(4);
Dcd = zeros(4,2);

% GenAugSysMatrix
[Psi, Gamma, Theta, Phi] = GenAugSysMatrix(Ad, Bd, Bwd, Czd, Np, Nu);
%%
[Psi_c, Gamma_c, Theta_c] = GenAugSysMatrix(Ad, Bd, Bwd, Ccd, Np, Nu);

% Hessian matrix
P = Theta' * QQ * Theta + RR;

% Initial condition
x0 = [0 0 0 0]';
u0 = [0 0]';
y0 = Cyd * x0;

%% Control Loop
Tend = length(ref)-Np;
QP_Opt = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');

% t = 0:Ts:((size(ref,1)-1)*Ts);
t = 0:Ts:(size(ref,1)*Ts);

x = x0;
u = u0;
y = y0;

X_REC = x;
Y_REC = y;
U_REC = u;

for k = 1 : length(t)-2
    disp(['k = ', num2str(k)]);
    y = Cyd * x;
    
    refSmpl = ref(k,:);
    refHrznCol = repmat(refSmpl, 1, Np)';
    
    y_prev = y;
    err = refHrznCol - (Psi * x + Gamma * u);
    q = - Theta' * QQ * err;
    
    [G, h] = GenConstMatrix(Bd, Ccd, Np, Nu, x, u, Psi_c, Gamma_c, Theta_c, u_min, u_max, z_min, z_max, du_min, du_max);
    P = (P + P')/2;
    du_QP = quadprog(P, q, G, h, [], [], [], [], [], QP_Opt);
    zz = Psi * x + Gamma * u + Theta * du_QP;

    z_Pred = reshape(zz, [nz, Np])';
    u_Pred = u + cumsum(u);
    
    du = reshape(du_QP, [nu, Nu])';
    
    z = zz(1:nz, :);
    u = u + du(1, :)';
    dis = disturbance(k,:)';

%     x = Ad * x + Bd * u;
%     x = Ad * x + Bd * (u + dis);
    x = Ad * x + Bd * u + Bd * dis;
    
    Y_REC = horzcat(Y_REC, y);
    U_REC = horzcat(U_REC, u);
    X_REC = horzcat(X_REC, x);
end

% Data Processing
T_REC = t(1:end-1);
X_REC(1, :) = X_REC(1, :) + h10_nmp;
X_REC(2, :) = X_REC(2, :) + h20_nmp;
X_REC(3, :) = X_REC(3, :) + h30_nmp;
X_REC(4, :) = X_REC(4, :) + h40_nmp;
U_REC(1, :) = U_REC(1, :) + v10_nmp;
U_REC(2, :) = U_REC(2, :) + v20_nmp;

REF_REC(1, :) = ref(:, 1)'/kc + h10_nmp;
REF_REC(2, :) = ref(:, 2)'/kc + h20_nmp;

%% Plotting
figure(1);
clf
subplot(3,2,1)
hold on
stairs(T_REC, X_REC(3,:),'b-', 'linewidth', 2)
ylabel('h_3 [cm]')
title('h3')
grid on;

subplot(3,2,2)
hold on
stairs(T_REC,X_REC(4,:), 'b-', 'linewidth', 2)
ylabel('h_4 [cm]')
title('h4')
grid on;

subplot(3,2,3)
hold on
stairs(T_REC,X_REC(1,:), 'b-', 'linewidth', 2)
stairs(T_REC, REF_REC(1, :), 'm-.', 'linewidth', 1)
ylabel('h_1 [cm]')
title('h1')
grid on;

subplot(3,2,4)
hold on
stairs(T_REC,X_REC(2,:), 'b-', 'linewidth', 2)
stairs(T_REC, REF_REC(2, :), 'm-.', 'linewidth', 1)
ylabel('h_2 [cm]')
title('h2')
grid on;

subplot(3,2,5)
hold on
stairs(T_REC,U_REC(1,:),' b-', 'linewidth', 2)
ylabel('u_1 [V]')
xlabel('t [s]')
grid on;

subplot(3,2,6)
hold on
stairs(T_REC,U_REC(2,:),' b-', 'linewidth', 2)
title('u2')
ylabel('u_2 [V]')
xlabel('t [s]')
zoom on
grid on;