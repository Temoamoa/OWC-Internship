%% Résolution de la pression dans la chambre de Mutriku, en incompressible

% Hypothèses : Fluide incompressible, turbine wells linéaire
% modèle de vague sans turbine dans l'owc (pas d'influence de la pression sur la
% colonne d'eau)

% On résout p(t) = (Ro*Omega/(K*D))*dV/dt
% Ref : On thermodynamics in the primary power conversion of oscillating
% water column wave energy converters
% Sheng, Wanan; Alcorn, Raymond; Lewis, Anthony

% dV/dt = débit volumique de l'air dans la chambre
% q = -dV/dt si Exhalation
% q = +dV/dt si Inhalation

function [Pression,PTO] = f_incompressible(t,SS,Vitesse_colonne)

load Constante.mat
load Wells.mat

S = Constante.S;
k1_turb = Constante.k1_turb(SS); %Omega = valeur pour chaque SS (sea state)
p_atm = Constante.p_atm;
Ro_atm = Constante.Ro_atm;

% Plot de quelques vagues dans l'OWC
% Débit volumique eau
dV_eau = Vitesse_colonne'*S; % vect ligne % Débit volumique eau OWC = dV/dt

%% Pression dans la chambre lorsque le gaz est incompressible - Turbine Wells Linéaire
p_incompressible = zeros(1,length(t));

for n = 1 : length(t)
    p_incompressible(n) = ((Ro_atm/p_atm)*k1_turb*dV_eau(n)); %adimmensionné
    %p_incompressible(n) = (Ro_atm*k1*dV_eau(n));    
end

[Pp, P_t, dm_t, eta] = PTO_Turb(p_incompressible,'inc',Constante,WellsParams,SS);

PTO = struct();
PTO.Pp = Pp;
PTO.P_t = P_t;
PTO.dm_t = dm_t;
PTO.eta = eta;

Pression = p_incompressible;

