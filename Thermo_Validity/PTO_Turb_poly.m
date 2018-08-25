function [Pp, P_t, dm_t, eta, valve, Psi, Phi, Tau] = PTO_Turb_poly(p_star,Omega,Turb,Constante)

% Constante
D = Turb.D;
p_atm = Constante.p_atm;
rho_atm = Constante.Ro_atm;
gamma = Constante.gamma;

%On le rajoutera après
%Ro_in
% if strcmp(Modele,'isen') || strcmp(Modele,'poly')
%     Ro_in = rho_atm*(max(p_star+1,1)).^(1/gamma);
% elseif strcmp(Modele,'inc')
%     Ro_in = rho_atm;
% end
rho_in = rho_atm*(max(p_star+1,1))^(1/gamma);
%rho_in = rho_atm;

%% Biradial
Omega = real(Omega);

% Pressure across the turbine
p_turbine = p_star * p_atm;

% Turbine characteristic curves
Psi = abs(p_turbine / (rho_in * D^2 * Omega^2 ));
        
if( Psi > Turb.Psi_max)
    Psi = Turb.Psi_max; 
end

% Flow rates
valve = 1;        
Phi = valve * sign(Psi) * 0.195 * (abs(Psi)+0.000001)^0.6;
   
Q_t = valve * sign( p_star ) * Omega *  Turb.D^3 * Phi;        
dm_t = rho_in * Q_t;       

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

% Turbine power and torque
P_t = eta *  p_atm * p_star * Q_t;
T_t = P_t / Omega; 
Tau=0;

Pp = p_turbine * Q_t;
end

