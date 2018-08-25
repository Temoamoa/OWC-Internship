%% Résolution de la pression dans la chambre de Mutriku

% Hypothèse : Air compressible, gaz parfait, isentropique
% turbine wells linéaire excitation
% modèle de vague dans l'owc avec wells turbine


function [Ro, p_gp_isen, PowerTurb] = f_gp_isen(t_simu,SS,y0,Position,Vitesse)

load Constante.mat
load Wells.mat

%Timing EDO -- Précision de l'edo
sim_dt = 0.05;
dt = sim_dt; % ON pourrait le réduire pour plus de précision 
t0 = t_simu(1);
tf_desire = t_simu(end);

t_edo = t0:dt:tf_desire;

% Résolution EDO
%Runge-Kutta ordre 4
y_RK4 = zeros(1,length(t_edo));
y_RK4(1) = y0; % IC

for n = 1:length(t_edo)-1

    %RK4
    k1 = Oceanet_Modl(t_edo(n),y_RK4(n),t_simu,Position,Vitesse,Constante,SS);
    k2 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k1*dt,t_simu,Position,Vitesse,Constante,SS);
    k3 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k2*dt,t_simu,Position,Vitesse,Constante,SS);
    k4 = Oceanet_Modl(t_edo(n)+dt,y_RK4(n)+k3*dt,t_simu,Position,Vitesse,Constante,SS);
    y_RK4(n+1) = y_RK4(n)+(dt/6)*(k1+2*k2+2*k3+k4); 
    
end

[Pp, P_t, dm_t, eta,Ro_in] = PTO_Turb(y_RK4,'isen',Constante,WellsParams,SS);

p_gp_isen = y_RK4';
Ro = Ro_in';

PowerTurb = struct();
PowerTurb.Pp = Pp;
PowerTurb.P_t = P_t;
PowerTurb.dm_t = dm_t;
PowerTurb.eta = eta;

end

% Modèle OCEANET, Gaz Parfait, isentropique
% Ref : OCEANET_Formulation et Polytropic_V1

function dpdt_star = Oceanet_Modl(t_edo,p_star,T_simu,Position,Vitesse,Constante,SS)

% Constante
S = Constante.S;
K = Constante.K;
D = Constante.D;
p_atm = Constante.p_atm;
Ro_atm = Constante.Ro_atm;
V0 = Constante.V0;
gamma = Constante.gamma;
Beta_g = Constante.Beta;
Omega = Constante.Omega(SS);

% Interpolation
position = interp1(T_simu',Position,t_edo);
vitesse = interp1(T_simu',Vitesse,t_edo);

% Volumes et débits
% Colonne d'eau
V_eau = position*S; % vect ligne % Volume OWC
dV_eau = vitesse*S; % vect ligne % Débit volumique eau OWC

V_air = V0 - V_eau;
dV_air = -dV_eau;

% PTO paramètre
Ro_in = Ro_atm*(max(p_star+1,1))^(1/gamma); 

% 1 - Dimensionless pressure head
psi = abs(p_star*p_atm/(Ro_in*D^2*Omega^2)) ; 

% 2 - Dimensionless flow rate
phi = K*psi; 

Qt = sign(p_star)*Omega*D^3*phi;
dm_t = Ro_in*Qt;

% Je pense que Ro_in est l'alternance entre Ro_chamber puis Ro_atm pour
% l'exhalation et l'inhalation ici dm_t = (K*D/Omega)*p_star*p_atm donc il
% ne se sert pas de Ro. Par conséquent Ro_in aura de l'influence que sur
% les paramètre de la turbine : psi,phi,eta et Pi

% Ma Fonction EDO 
dpdt_star = (-gamma *(p_star+1)*(dV_air/V_air)) ...
        - ((gamma*(p_star+1)^(Beta_g))*(dm_t/(Ro_atm*V_air))) ; 
end


