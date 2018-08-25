%% Init_Mutriku_Turb
%% Biradial turbine
% clc
% clear
% 
% n_SS = 5;
% a = Initiat_Biradial(n_SS);

function Turb = Init_Biradial_var(vect_Omega)

Turb = struct();
%n_SS = 10; %Sea state

Turb.Psi_max = 0;
Turb.SP = zeros(4,4);
% Turb.Psi_tbl = zeros(46,1);
% Turb.Phi_tbl = zeros(46,1);
% Turb.Eta_tbl = zeros(46,1);

Turb.D = 0.5;% Turbine diameter
Turb.I = 5;  % Turbine inertia 

% Turbine polynomial coefficients for efficiency calculation
Turb.Psi_max = 80;
Turb.SP = [-0.00311883, 0.0416233, -0.218747, 0.851294;
            15.8815, -17.9453, 6.53643, 0.00932473;
            28.2129, -29.0032, 9.83849, -0.319012;
            164.414, -101.898, 22.8393, -1.0917];

% Generator
Turb.Pnom = 30;                % Nom power [kW]
Pnom = Turb.Pnom *1000;        % Nom power [W]
Om_nom= 157.0796;             % Nom rot speed Om in [rad/s]
Nnom = Om_nom*30/pi;          % Nom rot speed N  in [rpm]
nb_pp = 2; % Number of pairs of pole

% Mechanical losses calculation
D_shaft = 2e-14*Pnom^3 - 7.e-9*Pnom^2 + 1e-3*Pnom+27.451; % Shaft diameter [m]
x1 = 6e-4*Nnom;
x2 = -1e-5*Nnom;
x3 = 2e-7*Nnom;
if Pnom < 300000
    Mech_loss = x3*D_shaft^3 + x2*D_shaft^2 + x1*D_shaft;
else
    Mech_loss = Pnom/300000 * (x3*237.451^3 + x2*237.451^2 + x1*237.451);
end

Turb.Pnom = Pnom;          % Nom power [W]
Turb.Om_nom= Om_nom;       % Nom rot speed Om in [rad/s]
Turb.Nnom = Nnom;          % Nom rot speed N  in [rpm]
Turb.Tnom = Pnom/Om_nom;   % Nom torque [Nm]
Turb.MmMn = 2;             % Ratio torque max/min
Turb.F_hz = Nnom*nb_pp/60; % Nominal frequency [Hz]
Turb.Vnom = 400;           % Nominal tension 400
Turb.Vmax = 690;           % Max voltage from power converter
Turb.over_speed = Turb.Vmax/Turb.Vnom; % Overspeed allowed by ratio of voltage from PE / generator
Turb.Vslope = Turb.Vnom/Turb.F_hz; % Slope of voltage actualisation for v = Vslope*f_hz
Turb.Tmax = Turb.Tnom * Turb.MmMn; % Max Torque [W]
Turb.Pmax = Turb.Tmax*Turb.Om_nom*Turb.over_speed; % Max power [W]
Turb.Om_max = 3000 *pi/30;% Maximal turbine speed 
Turb.Om_co =  2500 *pi/30;% Turbine cut-off speed (close valve) 
Turb.Om_ci =  1800 *pi/30;% Turbine cut-in  speed (open  valve) 
Turb.T_min = 0.1 *Turb.Tnom; % Electrical torque when valve is closed
Turb.nb_pp = nb_pp;        % Number of pairs of poles
Turb.Mech_loss = Mech_loss;% Mechanical losses [W]
Turb.Weight = 250;         % Weight [kg]
Turb.R_ohm = 72584*Pnom^(-1.236); % (Stator/rotor?) resistance [Ohm]
Turb.Kk = 44.8068;         % Eddy current loss coefficient
Turb.Kh = 5.0098;          % Hysteresis constant
Turb.Thick_max = 0.64;     % Core plate max thickness [m]
Turb.Flux_dens = 0.8;      % Flux density [Wb]      

Turb.Om_min = 50;
Turb.Om_Tci = 60;

%% Rotational speed control
% Fixed speed
Turb.vect_Omega = vect_Omega;%[110 110 120 140 160 180 200 210 240 260 280 290 310 310];
% for i = 1:14
%     if Turb.Om_vect(i) > Turb.Om_co
%         Turb.Om_vect(i) = Turb.Om_co*0.96;              
%     end
% end 
% Turb.Om_opt = Turb.Om_vect(n_SS);

% Variable speed
% Choose control algo
% 0: Batch BeP
% 1: FS PI-ctrl
% 2: VS TorqueLaw k_opt
% 3: VS TorqueLaw best fit
% 4: MPC
n_CL = 2;
switch n_CL
    case 2
        % BeP controller
%                 [eta, idx] = max(Turb.Eta_tbl);
%                 Psi_opt = Turb.Psi_tbl(idx);
%                 Phi_opt = Turb.Phi_tbl(idx);
%                 Kopt = eta * Psi_opt * Phi_opt * WEC.rho_at * Turb.D^5;
        Turb.a = 1.1e-3; Turb.b = 2;
    otherwise
        % Best fit controller
        Turb.a =   0.0009241;%  (7.603e-05, 0.001772)
        Turb.b =       2.046;% 
        % Turb.a = 0.000147265244267;
        % Turb.b = 2.37065478219;
end        

Turb.k_sigpsi =  607.378*1;        
Turb.Window_size = 180;
end