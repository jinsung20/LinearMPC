function [G, h] = GenConstMatrix(Bd, Ccd, Np, Nu, x0, u0, Psi_c, Gamma_c, Theta_c, u_min, u_max, z_min, z_max, du_min, du_max)
nu = size(Bd, 2);
nc = size(Ccd, 1);

G = zeros(2*Nu*nu + 2*Nu*nu + 2*Np*nc, Nu*nu);
h = zeros(2*Nu*nu + 2*Nu*nu + 2*Np*nc, 1);

G1 = zeros(2*Nu*nu, Nu*nu);
h1 = zeros(2*Nu*nu, 1);

%-------------------------------------
% for du constraints matrix-----------
%-------------------------------------
% for G matrix (LHS for du constraints)
for i = 1: Nu * nu
    G1(2*(i-1)+1 : 2*(i-1)+2, i) = [1; -1];
end

for i = 1 : Nu
    for ii = 1 : nu
        h1((2*nu*(i-1) + 1) + 2*(ii-1) : 2*nu*(i-1) + 2*ii) = [du_max(ii); -du_min(ii)];
    end
end

%-------------------------------------
% for u constraints matrix------------
%-------------------------------------
% for F0
F0 = zeros(2*Nu*nu, Nu*nu);
tF0 = kron(eye(nu), [1; -1]);

for i = 1 : Nu
%     F0(2*nu*(i-1)+1 : 2*nu*i, 1:nu) = tF0;
    F0(2*nu*(i-1)+1 : 2*nu*i, nu*(i-1)+1 : nu*(i-1)+nu) = tF0;
end

F1 = F0(:, 1:nu);

% for f
f = zeros(2*Nu*nu, 1);
for i = 1 : Nu
    for ii = 1 : nu
        f((2*nu*(i-1) + 1) + 2*(ii-1) : 2*nu*(i-1) + 2*ii) = [u_max(ii); -u_min(ii)];
    end
end

uConRhs = -F1 * u0 + f;
uConRhs = uConRhs(:, 1);

% G2 = F0
% h2 = uConRhs

%-------------------------------------
% for state constraints matrix--------
%-------------------------------------
StatePred = Psi_c * x0 + Gamma_c * u0;

% for G
% Initial Test Code
zG = kron(eye(Np*nc), [1; -1]);

zg = zeros(2*Np*nc, 1);
for i = 1 : Np
    for ii = 1 : nc
        zg((2*nc*(i-1) + 1) + 2*(ii-1) : 2*nc*(i-1) + 2*ii) = [z_max(ii); -z_min(ii)];
    end
end

zG_Final = zG * Theta_c;
zg_Final = zg - zG * StatePred;

% G = [G1; F0; zG_Final];
% h = [h1; uConRhs; zg_Final];

G(1:2*Nu*nu, 1:Nu*nu, :) = G1;
G(2*Nu*nu + 1 : 2*Nu*nu + 2*Nu*nu, :) = F0;
G(2*Nu*nu + 2*Nu*nu + 1 : 2*Nu*nu + 2*Nu*nu + 2*Np*nc, :) = zG_Final;

h(1:2*Nu*nu, :) = h1;
h(2*Nu*nu + 1 : 2*Nu*nu + 2*Nu*nu, :) = uConRhs;
h(2*Nu*nu + 2*Nu*nu + 1 : 2*Nu*nu + 2*Nu*nu + 2*Np*nc, :) = zg_Final;