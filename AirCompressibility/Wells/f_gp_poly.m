%% R�solution de la pression dans la chambre de Mutriku

% Hypoth�se : Air compressible, gaz parfait, isentropique
% turbine wells lin�aire excitation
% mod�le de vague dans l'owc avec wells turbine


function [Ro, p_gp_poly, PowerTurb] = f_gp_poly(t_simu,SS,y0,Position,Vitesse)

load Constante.mat
load Wells.mat

%Timing EDO -- Pr�cision de l'edo
sim_dt = 0.05;
dt = sim_dt*5; % ON pourrait le r�duire pour plus de pr�cision 
t0 = t_simu(1);
tf_desire = t_simu(end);

t_edo = t0:dt:tf_desire;

% R�solution EDO
%Runge-Kutta ordre 4
y_RK4 = zeros(1,length(t_edo));
y_RK4(1) = y0;% IC

for n = 1:length(t_edo)-1

    %RK4
    k1 = Oceanet_Modl(t_edo(n),y_RK4(n),t_simu,Position,Vitesse,Constante,WellsParams,SS);
    k2 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k1*dt,t_simu,Position,Vitesse,Constante,WellsParams,SS);
    k3 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k2*dt,t_simu,Position,Vitesse,Constante,WellsParams,SS);
    k4 = Oceanet_Modl(t_edo(n)+dt,y_RK4(n)+k3*dt,t_simu,Position,Vitesse,Constante,WellsParams,SS);
    y_RK4(n+1) = y_RK4(n)+(dt/6)*(k1+2*k2+2*k3+k4); 
    
end

[Pp, P_t, dm_t, eta,Ro_in] = PTO_Turb(y_RK4,'poly',Constante,WellsParams,SS);

p_gp_poly = y_RK4';
Ro = Ro_in';

PowerTurb = struct();
PowerTurb.Pp = Pp;
PowerTurb.P_t = P_t;
PowerTurb.dm_t = dm_t;
PowerTurb.eta = eta;

end

% Mod�le OCEANET, Gaz Parfait, Polytropique
% Ref : OCEANET_Formulation et Polytropic_V1

function dpdt_star = Oceanet_Modl(t_edo,p_star,T_simu,Position,Vitesse,Constante,WellsParams,SS)

% Constante
S = Constante.S;
Ro_atm = Constante.Ro_atm;
V0 = Constante.V0;
Beta_g = Constante.Beta;

% Interpolation
position = interp1(T_simu',Position,t_edo);
vitesse = interp1(T_simu',Vitesse,t_edo);

% Volumes et d�bits
% Colonne d'eau (Volume d'eau relatif
V_eau = position*S; % vect ligne % Volume OWC
dV_eau = vitesse*S; % vect ligne % D�bit volumique eau OWC

V_air = V0 - V_eau;
dV_air = -dV_eau;

% PTO param�tre 

[Pp, P_t,dm_t, eta] = PTO_Turb(p_star,'poly',Constante,WellsParams,SS);

% Kappa
%Isen
%Kappa = gamma;
%Beta_k = Beta_g;

%Poly
Kappa = 1/(1-(Beta_g*eta));
Beta_k = (Kappa-1)/Kappa;

% Je pense que Ro_in est l'alternance entre Ro_chamber puis Ro_atm pour
% l'exhalation et l'inhalation ici dm_t = (K*D/Omega)*p_star*p_atm donc il
% ne se sert pas de Ro. Par cons�quent Ro_in aura de l'influence que sur
% les param�tre de la turbine : psi,phi,eta et Pi

% Ma Fonction EDO 
dpdt_star = (-Kappa *(p_star+1)*(dV_air/V_air)) ...
        - ((Kappa*(p_star+1)^(Beta_k))*(dm_t/(Ro_atm*V_air))) ; 
end


