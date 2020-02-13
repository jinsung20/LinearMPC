function [Psi, Gamma, Theta, Phi] = GenAugSysMatrix(A, B, Bw, Cz, Np, Nu)  

nx = size(A, 2);
nu = size(B, 2);
nw = size(Bw, 2);
nz = size(Cz, 1);

% -------------------------------------
%  for Psi ----------------------------
% -------------------------------------
% CAterm = Cz;
% Psi = zeros(0, size(A,2));
% for i = 1 : Np
%     CAterm = CAterm * A;
%     Psi = vertcat(Psi, CAterm);
% end
% Psi_orig = Psi;


CAterm = Cz;
Psi = zeros(nz * Np, nx);
for i = 1 : Np
    CAterm = CAterm * A;
    Psi((i-1)*nz + 1 : i*nz, :) = CAterm;
end


% -------------------------------------
%  for Gamma --------------------------
% -------------------------------------
Gamma = zeros(nz * Np, nu);
tGamma = zeros(nz, nu);
for i = 1 : Np
    Gamma((1+(i-1)*nz):i*nz, :) = tGamma + Cz*A^(i-1)*B;
    tGamma = Gamma((1+(i-1)*nz):i*nz, :);
end


% -------------------------------------
%  for Theta --------------------------
% -------------------------------------
Theta = zeros(nz * Np, nu * Nu);
Theta(:, 1:nu) = Gamma;

t1 = zeros(nz * Np, nu);
for i = 1 : Nu - 1
    t1(1:i*nz, :) = zeros(size(t1(1:i*nz, :)));
    t1(i*nz+1:end, :) = Theta(1:size(Theta,1)-(i*nz), 1:nu);
    Theta(:, i*nu+1:(i+1)*nu) = t1;
end


% -------------------------------------
%  for Phi --------------------------
% -------------------------------------
Phi = zeros(Np*nz, nw);
tPhi = zeros(nz, nw);
for i = 1 : Np
    Phi((1+(i-1)*nz):i*nz, :) = tPhi + Cz*A^(i-1)*Bw;
    tPhi = Phi((1+(i-1)*nz):i*nz, :);
end
