clc
clearvars
close all

%load Data\Constante.mat
load BD_Birad_SS.mat
load Constante.mat

Turbine_type = 3;
CL = 1;

Position_birad = BD_Birad_SS.Position_birad;
Vitesse_birad = BD_Birad_SS.Vitesse_birad;
Omega = BD_Birad_SS.Omega;
p_star = BD_Birad_SS.p_star;

V_eau = Position_birad*Constante.S; % vect ligne % Volume OWC
dV_eau = Vitesse_birad*Constante.S; % vect ligne % Débit volumique eau OWC

V_air = -V_eau; %Relative volume %Constante.V0 - V_eau; %Real value
dV_air = -dV_eau;

t0 = 0;
tf_desire = 1800;
sim_dt = 0.05 ; %incrément de temps
%t = t0:sim_dt:tf_desire ;
t = linspace(t0,tf_desire,length(Omega));
t_edo = t0:sim_dt*5:tf_desire ; % a multiplié ou divisé sim_dt pour accéléré ou avoir plus de précision
%t_edo = t0:sim_dt*5:tf_desire ;

if t(1) == 0
    t_ini1 = 1;
    t_end1 = length(t);   
else
    t_ini1 = t(1)/sim_dt;
    t_end1 = t_ini1 + length(t) - 1;
end

Indice_temp = (t_ini1 : t_end1)';


Y0 = p_star(Indice_temp(1));

%% Décalage entre flow rate et pression

% for SS = 13%1:14
%     [WEC, Turb] = Init_Mutriku(Turbine_type,SS,CL);
%     [p_poly] = f_gp_poly_birad(t,Y0,t_edo,Position_birad(Indice_temp,SS),Vitesse_birad(Indice_temp,SS),Omega(Indice_temp,SS),WEC,Turb);
%     plot(t_edo,p_poly,t,dV_eau(Indice_temp,SS)/100)
%     legend('Pressure','Flow rate')
%     title(sprintf('Biradial, Polytropic model, CL1, SS = %d',SS))
%     xlabel('t[s]')
%     ylabel('q/100 [m^3.s^{-1}], p* [-]')
% end

%% Err relative entre Poly et Isen au niveau du Birad

ErrBirad_PolyIsen2 = zeros(14,1);
%figure()
for SS = 1%:14
    [WEC, Turb] = Init_Mutriku(Turbine_type,SS,CL);
    [p_poly] = f_gp_poly_birad(t,Y0,t_edo,Position_birad(Indice_temp,SS),Vitesse_birad(Indice_temp,SS),Omega(Indice_temp,SS),WEC,Turb); %,PowerTurb_poly
    [p_isen] = f_gp_isen_birad(t,Y0,t_edo,Position_birad(Indice_temp,SS),Vitesse_birad(Indice_temp,SS),Omega(Indice_temp,SS),WEC,Turb);
    
    ErrBirad_PolyIsen2(SS) = 100 * ( sum(abs(p_isen)) - sum(abs(p_poly)) ) ./ sum(abs(p_isen)) ;
    
    %Plot Intéractif des pression poly VS isen VS inc
    %hold on
%     plot(t_edo,p_poly,t_edo,p_isen)%,t_edo,zeros(1,length(t_edo)))
%     legend('poly','isen')
%     title(sprintf('Biradial, Sea state = %d',SS))
%     xlabel('t[s]')
%     ylabel('p*[-]')
%     grid on
%     drawnow
end
%%
%save 'Err_Birad2.mat' ErrBirad_PolyIsen2