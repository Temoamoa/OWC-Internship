clc
clearvars
close all

%load Data\Constante.mat
load BD_simulink_OWC.mat
load Constante.mat

% load Data\DBReality_biradial.mat
% Position_birad = DB_Biradial.DB{7,1}.Position_air*Constante.S;
% Position_birad = Position_birad - DB_Biradial.DB{7,1}.H0*Constante.S; % Relative Position
% Vitesse_birad = DB_Biradial.DB{7,1}.Vitesse_air'*Constante.S;
% p_data = DB_Biradial.DB{7,1}.pPort';
% p_star_data  = p_data/Constante.p_atm;

Position_wells = BD_simulink_OWC.Position_wells;
Vitesse_wells = BD_simulink_OWC.Vitesse_wells;
p_star_FX = BD_simulink_OWC.p_star_FX;

V_eau = Position_wells*Constante.S; % vect ligne % Volume OWC
dV_eau = Vitesse_wells*Constante.S; % vect ligne % Débit volumique eau OWC

V_air = -V_eau; %Relative volume %Constante.V0 - V_eau; %Real value
dV_air = -dV_eau;


t0 = 100;
tf_desire = 200;
sim_dt = 0.05 ; %incrément de temps
t = t0:sim_dt:tf_desire ;
t_edo = t0:sim_dt*5:tf_desire ;

if t(1) == 0
    t_ini1 = 1;
    t_end1 = length(t);   
else
    t_ini1 = t(1)/sim_dt;
    t_end1 = t_ini1 + length(t) - 1;
end

Indice_temp = (t_ini1 : t_end1)';


%% Perte de Puissance Pneumatique et Turbine

Y0 = p_star_FX(Indice_temp(1));
Pneuma = zeros(1,14);
Turb = zeros(1,14);
haut_isen = zeros(1,14); 
haut_poly = zeros(1,14);
bas_isen = zeros(1,14);
bas_poly = zeros(1,14);

for SS = 13%1:14
    sprintf('Sea state = %d',SS)
    [p_inc,PowerTurb_inc] = f_incompressible(t,SS,Vitesse_wells(Indice_temp,SS));
    figure()
    subplot(2,1,1)
    plot(t,p_inc,t,dV_eau(Indice_temp,SS)/100)
    legend('Pressure','Flow rate')
    title(sprintf('Wells, Incompressible model, SS = %d',SS))
    xlabel('t[s]')
    ylabel('q/100 [m^3.s^{-1}], p* [-]')
    grid on 
    
    [Ro_poly,p_poly,PowerTurb_poly] = f_gp_poly(t,SS,Y0,Position_wells(Indice_temp,SS),Vitesse_wells(Indice_temp,SS));
    subplot(2,1,2)
    plot(t_edo,p_poly,t,dV_eau(Indice_temp,SS)/100)
    legend('Pressure','Flow rate')
    title(sprintf('Wells, Polytropic model, SS = %d',SS))
    xlabel('t[s]')
    ylabel('q/100 [m^3.s^{-1}], p* [-]')
    grid on
end

%%
figure()
for SS = 13%1:14
    [p_inc,PowerTurb_inc] = f_incompressible(t,SS,Vitesse_wells(Indice_temp,SS));
    [Ro_poly,p_poly,PowerTurb_poly] = f_gp_poly(t,SS,Y0,Position_wells(Indice_temp,SS),Vitesse_wells(Indice_temp,SS));
    [Ro_isen,p_isen,PowerTurb_isen] = f_gp_isen(t,SS,Y0,Position_wells(Indice_temp,SS),Vitesse_wells(Indice_temp,SS));

    % PTO
    %Inc
    Pp_inc = PowerTurb_inc.Pp;
    P_t_inc = PowerTurb_inc.P_t;
    dm_t_inc = PowerTurb_inc.dm_t;
    eta_inc = PowerTurb_inc.eta;

    %Isen
    Pp_isen = PowerTurb_isen.Pp;
    P_t_isen = PowerTurb_isen.P_t;
    dm_t_isen = PowerTurb_isen.dm_t;
    eta_isen = PowerTurb_isen.eta;

    %Poly
    Pp = PowerTurb_poly.Pp;
    P_t = PowerTurb_poly.P_t;
    dm_t = PowerTurb_poly.dm_t;
    eta = PowerTurb_poly.eta;
    
    %Perte puissance comparaison entre modèle
    a = sum(Pp_inc(2:end));
    b = sum(Pp_isen(2:end));
    c = sum(Pp(2:end));
    Perte_Puissance_Pneumatique_isen = 100*abs((a-b)/a);
    Perte_Puissance_Pneumatique_poly = 100*abs((a-c)/a);
    Pneuma(SS) = abs(Perte_Puissance_Pneumatique_isen - Perte_Puissance_Pneumatique_poly);

    a = sum(P_t_inc(2:end));
    b = sum(P_t_isen(2:end));
    c = sum(P_t(2:end));
    Perte_Puissance_Turbine_isen = 100*abs((a-b)/a);
    Perte_Puissance_Turbine_poly = 100*abs((a-c)/a);
    Turb(SS) = abs(Perte_Puissance_Turbine_isen - Perte_Puissance_Turbine_poly);
    
    %Perte de pression partie haute
    a = max(p_inc);
    b = max(p_isen);
    c = max(p_poly);

    haut_isen(SS) = 100*abs((a-b)/a);
    haut_poly(SS) = 100*abs((a-c)/a);
    
    %Perte de pression partie basse
    a = min(p_inc);
    b = min(p_isen);
    c = min(p_poly);

    bas_isen(SS) = 100*abs((a-b)/a);
    bas_poly(SS) = 100*abs((a-c)/a);
    
    %Plot Intéractif des pression poly VS isen VS inc
    %hold on
    plot(t_edo,p_poly,t,p_isen,t,p_inc)%,t,zeros(1,length(t)))
    legend('poly','isen','inc')
    title(sprintf('Wells, Sea state = %d',SS))
    xlabel('t[s]')
    ylabel('p*[-]')
    grid on
    drawnow
end

%%
% Erreur puissance entre isen et poly
Puissance_perdu = [(1:14)' Pneuma' Turb'];

% Diff des modèles pressions
ERR = [(1:14)' haut_isen' bas_isen'  haut_poly' bas_poly'];
% Diff Isen-inc et poly-inc
ERRT = [(1:14)' (haut_isen'+bas_isen') (haut_poly'+bas_poly')];
% Diff Isen-Poly
ERR_Isen_Poly = [(1:14)' abs((haut_isen'+bas_isen')-(haut_poly'+bas_poly'))];

%%
PerteNRJtf1799 = struct();
PerteNRJtf1799.Puissance_perdu = Puissance_perdu;
PerteNRJtf1799.ERR = ERR;
PerteNRJtf1799.ERRT = ERRT;
PerteNRJtf1799.ERR_Isen_Poly = ERR_Isen_Poly;

%save('perteNRJtf1790.mat')

%% On peut constater que lorsque l'air est considéré incompressible il y a plus d'nrj
% %%
% SS = 10;
% [p_inc,PowerTurb_inc] = f_incompressible(t,SS,Vitesse_wells(Indice_temp,SS));
% %%
% [Ro_poly,p_poly,PowerTurb_poly] = f_gp_poly(t,SS,Position_wells(Indice_temp,SS),Vitesse_wells(Indice_temp,SS));
% [Ro_isen,p_isen,PowerTurb_isen] = f_gp_isen(t,SS,Position_wells(Indice_temp,SS),Vitesse_wells(Indice_temp,SS));
% 
% %% PTO
% %Inc
% Pp_inc = PowerTurb_inc.Pp;
% P_t_inc = PowerTurb_inc.P_t;
% dm_t_inc = PowerTurb_inc.dm_t;
% eta_inc = PowerTurb_inc.eta;
% 
% %Isen
% Pp_isen = PowerTurb_isen.Pp;
% P_t_isen = PowerTurb_isen.P_t;
% dm_t_isen = PowerTurb_isen.dm_t;
% eta_isen = PowerTurb_isen.eta;
% 
% %Poly
% Pp = PowerTurb_poly.Pp;
% P_t = PowerTurb_poly.P_t;
% dm_t = PowerTurb_poly.dm_t;
% eta = PowerTurb_poly.eta;
% 
% l = 2;
% c = 2;
% n = 1;
% 
% figure()
% subplot(l,c,n)
% plot(t,dm_t_inc,t,dm_t)
% legend('inc','poly')
% title('dmdt')
% 
% n = n+1;
% subplot(l,c,n)
% plot(t,Pp_inc,t,Pp)
% legend('inc','poly')
% title('Pp')
% 
% n = n+1;
% subplot(l,c,n)
% plot(t,P_t_inc,t,P_t)
% legend('inc','poly')
% title('Pt')
% 
% n = n+1;
% subplot(l,c,n)
% plot(t,eta_inc,t,eta)
% legend('inc','poly')
% title('eta')
% 
% %% Puissance pneumatique Vs Turbine pour le poly
% figure()
% plot(t,Pp,t,P_t)
% legend('Pneuma','Turb')
% title('Puissance pneumatique Vs Turbine pour le poly')

