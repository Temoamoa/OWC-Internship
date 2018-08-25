function [Pp, P_t, dm_t, eta,Ro_in] = PTO_Birad2(p_star,Omega,Modele,WEC,Turb)

% Constante

D = Turb.D;
p_atm = WEC.p_at;
Ro_atm = WEC.rho_at;
gamma = WEC.gamma;
Beta_g = WEC.beta;

%Ro_in
if strcmp(Modele,'isen') || strcmp(Modele,'poly')
    Ro_in = Ro_atm*(max(p_star+1,1)).^(1/gamma);
end

% Turbine characteristic curves of :
% 1 - Dimensionless pressure head
Psi = abs(p_star*p_atm./(Ro_in*D^2*Omega^2)) ; 

% 2 - Dimensionless flow rate
phi = sign(Psi)*0.195 * (abs(Psi)+0.000001)^0.6;

% 3 - Efficiency
eta = efficiency(Psi,Turb);

if strcmp(Modele,'poly')
    %Poly
    Kappa = 1./(1-(Beta_g*eta));
    %Correction Ro_in, phi, eta
    Ro_in = Ro_atm*(max(p_star+1,1)).^(1./Kappa); 
    Psi = abs(p_star*p_atm./(Ro_in*D^2*Omega^2)) ; 
    phi = sign(Psi)*0.195 * (abs(Psi)+0.000001)^0.6;
    eta = efficiency(Psi,Turb);
end

% Flow rates with dimensions
Qt = sign(p_star)*Omega*D^3*phi;
dm_t = Ro_in.*Qt;

Pp = dm_t.*p_star*p_atm./Ro_in;

% Turbine power and torque
T_t = Ro_in.*Omega^2.*D^5.*eta.*phi.*Psi;
P_t = T_t*Omega;

% The dimensionless power coefficient
%Pi = P_t/(Ro_in*Omega^3*D^5);
end

function eta = efficiency(Psi,Turb)
if( Psi > 3.3116337496519517 )
        eta = 1.0/(1.442750893460873 + 0.2066910461657107*Psi); 
    elseif( Psi > 0.35075016168255563 )
        idx = 1;
        eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
    elseif( Psi > 0.3082881279543439 )
        idx = 2;
        eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
    elseif( Psi > 0.18131721834023762 )
        idx = 3;
        eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
    else
        idx = 4;
        eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
end 

if isnan(eta), eta=0; end
end