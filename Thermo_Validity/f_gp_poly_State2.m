%% R�solution de la pression dans la chambre de Mutriku
%% Biradial

function [p_gp_poly] = f_gp_poly_State2(Y0,t_simu,dt_in,State,Vect_Omega,Position,Vitesse,V0)

load Constante.mat

% Initialisation Biradial
Turb = Init_Biradial_var(Vect_Omega);

%Simu data
sim_dt = dt_in; % de l'ordre de 0.25;

% Timing EDO -- Pr�cision de l'edo
dt = sim_dt; % ON pourrait le r�duire pour plus de pr�cision ou le multipli� pour que �a aille plus vite
t0 = t_simu(1);
tf_desire = t_simu(end);

t_edo = t0:dt:tf_desire;

% R�solution EDO
% IC
y0 = Y0;

%Runge-Kutta ordre 4
y_RK4 = zeros(1,length(t_edo));
y_RK4(1) = y0;

n = 1;
while n < length(t_edo)-1
    
    %state = interp1(t_simu,State',t_edo(n));
    
    if State(n) == 4
        k1 = Oceanet_Modl(t_edo(n),y_RK4(n),t_simu,Vect_Omega,Position,Vitesse,Turb,Constante);
        k2 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k1*dt,t_simu,Vect_Omega,Position,Vitesse,Turb,Constante);
        k3 = Oceanet_Modl(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k2*dt,t_simu,Vect_Omega,Position,Vitesse,Turb,Constante);
        k4 = Oceanet_Modl(t_edo(n)+dt,y_RK4(n)+k3*dt,t_simu,Vect_Omega,Position,Vitesse,Turb,Constante);
        y_RK4(n+1) = y_RK4(n)+(dt/6)*(k1+2*k2+2*k3+k4); 
        
    else %State(n) == 2
%         k1 = NoPTO(t_edo(n),y_RK4(n),t_simu,Position,Vitesse,Constante);
%         k2 = NoPTO(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k1*dt,t_simu,Position,Vitesse,Constante);
%         k3 = NoPTO(t_edo(n)+0.5*dt,y_RK4(n)+0.5*k2*dt,t_simu,Position,Vitesse,Constante);
%         k4 = NoPTO(t_edo(n)+dt,y_RK4(n)+k3*dt,t_simu,Position,Vitesse,Constante);
%         y_RK4(n+1) = y_RK4(n)+(dt/6)*(k1+2*k2+2*k3+k4);
        
        v0 = interp1(t_simu,V0',t_edo(n));
        volume = interp1(t_simu,Position',t_edo(n));
        y_RK4(n+1) = ((v0/volume)^Constante.gamma)-1;
    end

    n = n+1;
    
end

%[Pp, P_t, dm_t, eta, valve, Psi, Phi, Tau] = PTO_Turb(y_RK4,Turb,Constante);

p_gp_poly = y_RK4';
%Ro = Ro_in';

% PowerTurb = struct();
% PowerTurb.Pp = Pp;
% PowerTurb.P_t = P_t;
% PowerTurb.dm_t = dm_t;
% PowerTurb.eta = eta;
% PowerTurb.valve = valve;
% PowerTurb.Psi = Psi;
% PowerTurb.Phi = Phi;
% PowerTurb.Tau = Tau;
end

% Mod�le OCEANET, Gaz Parfait, Polytropique
% Ref : OCEANET_Formulation et Polytropic_V1

function dpdt_star = Oceanet_Modl(t_edo,p_star,T_simu,Vect_Omega,Position,Vitesse,Turb,Constante)
    % Interpolation
%     T0 = 293;
%     R = 8.31;
    
position = interp1(T_simu,Position',t_edo);
vitesse = interp1(T_simu,Vitesse',t_edo);
Omega = interp1(T_simu,Vect_Omega',t_edo);

V_air = position;
dV_air = vitesse;
    
%if State == 4
    % PTO param�tre 

[Pp, P_t, dm_t, eta] = PTO_Turb_poly(p_star,Omega,Turb,Constante);%, valve, Psi, Phi, Tau

%Poly
Kappa = 1/(1-(Constante.Beta*eta));
Beta_k = (Kappa-1)/Kappa;

% Je pense que Ro_in est l'alternance entre Ro_chamber puis Ro_atm pour
% l'exhalation et l'inhalation ici dm_t = (K*D/Omega)*p_star*p_atm donc il
% ne se sert pas de Ro. Par cons�quent Ro_in aura de l'influence que sur
% les param�tre de la turbine : psi,phi,eta et Pi

% Ma Fonction EDO 
dpdt_star = (-Kappa *(p_star+1)*(dV_air/V_air)) ...
        - ((Kappa*(p_star+1)^(Beta_k))*(dm_t/(Constante.Ro_atm*V_air))); 
    
% else
%     dpdt_star = 0;%(-Constante.gamma *(p_star+1)*(dV_air/V_air));
    %G = Constante.Ro_atm*dV_air;
    %T = @(p_star)(T0*(1+p_star/Constante.p_atm)^(Constante.Beta)-T0);
    %dpdt_star = ((R*G*(T0*(1+p_star/Constante.p_atm)^(Constante.Beta)-T0) - dV_air*p_star)/V_air) * (1/(1-(p_star*T0*(Constante.gamma-1))/((T0*(1+p_star/Constante.p_atm)^(Constante.Beta)-T0)*Constante.p_atm*Constante.gamma))*(1-(p_star/Constante.p_atm))^(-1/Constante.gamma));
end

function dpdt_star = NoPTO(t_edo,p_star,T_simu,Position,Vitesse,Constante)

position = interp1(T_simu,Position',t_edo);
vitesse = interp1(T_simu,Vitesse',t_edo);

V_air = position;
dV_air = vitesse;

dpdt_star = -Constante.gamma *(p_star+1)*(dV_air/V_air);

end