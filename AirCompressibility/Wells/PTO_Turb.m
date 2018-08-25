function [Pp, P_t, dm_t, eta,Ro_in] = PTO_Turb(p_star,Modele,Constante,WellsParams,SS)

% Constante
K = Constante.K;
D = Constante.D;
p_atm = Constante.p_atm;
Ro_atm = Constante.Ro_atm;
Omega = Constante.Omega(SS);
gamma = Constante.gamma;
Beta_g = Constante.Beta;

%Ro_in
if strcmp(Modele,'isen') || strcmp(Modele,'poly')
    Ro_in = Ro_atm*(max(p_star+1,1)).^(1/gamma);
elseif strcmp(Modele,'inc')
    Ro_in = Ro_atm;
end

% Turbine characteristic curves of :
% 1 - Dimensionless pressure head
psi = abs(p_star*p_atm./(Ro_in*D^2*Omega^2)) ; 

% 2 - Dimensionless flow rate
phi = K*psi; 

% 3 - Efficiency
Phi_tbl = WellsParams(:,1);
Phi_tbl(1) = 1e-10;
Eta_tbl = WellsParams(:,4);

% Turbine efficiency
eta = interp1(Phi_tbl,Eta_tbl,phi,'linear'); %linear is default
if isnan(eta), eta=0; end


if strcmp(Modele,'poly')
    %Poly
    Kappa = 1./(1-(Beta_g*eta));
    %Correction Ro_in, phi, eta
    Ro_in = Ro_atm*(max(p_star+1,1)).^(1./Kappa); 
    psi = abs(p_star*p_atm./(Ro_in*D^2*Omega^2)) ; 
    phi = K*psi; 
    eta = interp1(Phi_tbl,Eta_tbl,phi,'linear'); %linear is default
    if isnan(eta), eta=0; end
end

% Flow rates with dimensions
Qt = sign(p_star).*Omega.*D^3.*phi;
dm_t = Ro_in.*Qt;

Pp = dm_t.*p_star*p_atm./Ro_in;

% Turbine power and torque
T_t = Ro_in.*Omega^2.*D^5.*eta.*phi.*psi;
P_t = T_t*Omega;

% The dimensionless power coefficient
%Pi = P_t/(Ro_in*Omega^3*D^5);

end