%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Numercial modelling of Mutriku OWC                      %
%                                           OPERA Project                 %
%  F.X. Faÿ                                                               %
%  Jan 2017                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WEC, Turb] = Init_Mutriku_Data(DB_reality,Turbine_type,n_CL,mean_H0) %[WEC, Turb] = Init_Mutriku_Data(DB_reality,Turbine_type, n_CL)

% Physical parameters
rho_w = 1025;   % Water density
rho_at= 1.25;   % Air density
g     = 9.81;   % Gravitational constant
p_at  = 101500; % Atmospheric Pressure (Pa)
gamma = 1.4;    % Air Adiabatic expansion coefficient
beta  = 0.02857;% 

% Device parameters
w_ch = 4.50;              % Chamber width
l_ch = 3.10;              % Chamber length
h_ch = mean_H0; %7.45;              % Chamber height between B.M.V.E h_ch = 9.70 - 0.00 
                          %                        P.M.V.E h_ch = 9.70 - 4.50                         
S_iws = w_ch * l_ch;      % Area of the internal water surface (iws)
Vo = h_ch * S_iws ;       % Initial chamber air volume
K33= rho_w * g * S_iws;   % Hydrostatic stifness coefficient
Kp = p_at * S_iws;        % Pressure coefficient

%% Time domain analysis
% Hydrodynamic coefficients in FD
load Rad_SS.mat
order = length(Ar);

%% Matrices
% Equation of motion 
% (M+mu)*z_dotdot(t) = Fex(t) - Cr.Ir(t) + K33*z(t) + Kp * p*
% X_dot = F(t,Ir)
% X = [z_dot z p* Om Ir]' | U = Fex
%       X_dot = A * X + B * U + f(s)
%           Y = C * X + D * U
n_st = 4; % Number of states of the state vector (w/out rad) X = [z_dot, z, p*, Om]'
n_in = 1; % Number of inputs U = Fex
Zrx1 = zeros(order,1);
Z1xr = zeros(1,order);

M = (m+Mu);%*.5;

A = [ 0,  -K33/M,  -Kp/M,    0,  -Cr/M;
      1,       0,     0,    0,   Z1xr;
      0,       0,     0,    0,   Z1xr;
      0,       0,     0,    0,   Z1xr;
     Br,    Zrx1,  Zrx1, Zrx1,    Ar];

B = [1/M; zeros(n_st+order-1,n_in)];

C = [eye(4), zeros(n_st, order);
        zeros(order, order+4) ];

D = zeros(n_st+order,n_in);

%% Output data
% OWC Plant parameters
WEC.rho_w = rho_w; 
WEC.rho_at= rho_at; 
WEC.p_at  = p_at;
WEC.gamma = gamma;
WEC.beta  = beta;
WEC.h_ch  = h_ch;
WEC.S_iws = S_iws;
% Excitation force
WEC.Fint = Fint;
% State space matrices
WEC.A = A;
WEC.B = B;
WEC.C = C;
WEC.D = D;
% RAD matrices
WEC.Ar = Ar;
WEC.Br = Br;
WEC.Cr = Cr;
WEC.Dr = Dr;
WEC.order = order;
% Constants
WEC.K33 = K33;
WEC.Kp  = Kp;
WEC.Mu  = Mu;

% Initial condition
WEC.Omega0 = 100;
WEC.X0 = complex(zeros(24,1) ,0);%zeros(n_st+order,1);
WEC.X0(4) = complex(WEC.Omega0 ,0);

%% Turbine characteristics
Turb.Turb = Turbine_type;
% initialise specific parameters
Turb.Psi_max = 0;
Turb.SP = zeros(4,4);
Turb.Psi_tbl = zeros(46,1);
Turb.Phi_tbl = zeros(46,1);
Turb.Eta_tbl = zeros(46,1);

switch Turbine_type    
    case 1
        %% Wells Mutriku
        Gen_type = 1; % 1: 18.5kW 460V | 30kW 690V
        load Wells.mat
        Turb.D = 0.75;
        Turb.I = 3.06;
        
%   plot(WellsParams(:,1),WellsParams(:,3));hold on; plot(WellsParams(:,1),WellsParams(:,2));plot(WellsParams(:,1),WellsParams(:,4));hold off
        Turb.Kt = mean(WellsParams(4:end,1)./WellsParams(4:end,2));
        Turb.Phi_tbl = WellsParams(:,1);
        Turb.Tau_tbl = WellsParams(:,3);
        Turb.Eta_tbl = WellsParams(:,4); 
       
        % Generator characteristics
        switch Gen_type
            case 1
                Pnom = 18500; 
                Nnom = 3000;
                Turb.MmMn = 2;             % Ratio torque max/min
                Turb.Vnom = 460;           % Nominal tension 400
                Turb.Vmax = 460;           % Max voltage from power converter                
                Turb.Weight = 100;         % Weight [kg]
                Turb.Om_max = 4800 *pi/30;% Maximal turbine speed 
                Turb.Om_co =  4380 *pi/30;% Turbine cut-off speed (close valve) 
                Turb.Om_os =  3900 *pi/30;% Overspeed
                Turb.Om_ci =  3500 *pi/30;% Turbine cut-in  speed (open  valve) 80% Om_co
            case 2
                Pnom = 30000;
                Nnom = 1500;
                Turb.MmMn = 2;             % Ratio torque max/min
                Turb.Vnom = 400;           % Nominal tension 400
                Turb.Vmax = 690;           % Max voltage from power converter
                Turb.Weight = 280;         % Weight [kg]
                Turb.Om_max = 3000 *pi/30;% Maximal turbine speed 
                Turb.Om_co =  2500 *pi/30;% Turbine cut-off speed (close valve) 
                Turb.Om_ci =  1800 *pi/30;% Turbine cut-in  speed (open  valve) 
        end
        Turb.Om_min = 40;
        Turb.Om_Tci = 60;
          
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
        
        Turb.Pnom = Pnom;
        Turb.Nnom = Nnom;
        Turb.Om_nom = Turb.Nnom * pi/30;
        Turb.Tnom = Turb.Pnom/Turb.Om_nom;   % Nom torque [Nm]
        Turb.F_hz = Turb.Nnom/60; % Nominal frequency [Hz]   
        Turb.over_speed = Turb.Vmax/Turb.Vnom; % Overspeed allowed by ratio of voltage from PE / generator
        Turb.Vslope = Turb.Vnom/Turb.F_hz; % Slope of voltage actualisation for v = Vslope*f_hz
        Turb.Tmax = Turb.Tnom * Turb.MmMn; % Max Torque [W]
        Turb.Pmax = Turb.Tmax*Turb.Om_nom*Turb.over_speed; % Max power [W]
        Turb.T_min = 0.1 *Turb.Tnom; % Electrical torque when valve is closed
        Turb.nb_pp = 1;        % Number of pairs of poles
        Turb.Mech_loss = Mech_loss;% Mechanical losses [W]
        
        Turb.R_ohm = 72584*Turb.Pnom^(-1.236); % (Stator/rotor?) resistance [Ohm]
        Turb.Kk = 44.8068;         % Eddy current loss coefficient
        Turb.Kh = 5.0098;          % Hysteresis constant
        Turb.Thick_max = 0.64;     % Core plate max thickness [m]
        Turb.Flux_dens = 0.8;      % Flux density [Wb]
                
        Turb.Psi_max = 0;
        Turb.SP = zeros(4,4);

        %% Rotational speed control
        % Fixed speed
        Turb.Om_vect = [80 95 100 100 110 120 130 140 155 160 175 180 210 215];%[100 100 110 130 150 170 190 210 240 250 280 290 310 320];%[80 80 80 90 90 100 110 110 120 130 140 140 160 160];%[100 100 110 125 140 165 180 195 220 230 255 265 280 290];
        for i = 1:14
            if Turb.Om_vect(i) > Turb.Om_co
                Turb.Om_vect(i) = Turb.Om_co*0.96;              
            end
        end
        Turb.Om_opt = Turb.Om_vect(n_SS);
                
        % Variable speed
        switch n_CL
            case 2
                % BeP controller
                Kopt_wells = [0.000170806000000000 4.73970000000000e-05 0.000310828100000000 0.00546579140000000 0.00151670350000000 0.00994650070000000];
                Tau_max = max(Turb.Tau_tbl);
                Kopt = Tau_max * WEC.rho_at * (Turb.D/2)^5;
%                 [Eta_max idx] = max(Turb.Eta_tbl);
%                 Kopt = Eta_max * Turb.Phi_tbl(idx)*(1+1/Turb.Kt) * WEC.rho_at * (Turb.D/2)^5;
%                 Kopt = 3.88655E-05;%with D/2 | 0.0012436969;with D % Phi = 0.0785861	Psi = 0.147872  Eta = 0.44 eta_avg in CL1 left of curve
                Turb.a = Kopt; Turb.b = 2;           
%                 Turb.a = Kopt_wells(6);

            case 3
                % Best fit controller
                Turb.a =    3.348e-06;%  Opt 032018
                Turb.b =       3.099 ;

            otherwise
                Turb.a =    3.348e-06;%  Opt 032018
                Turb.b =       3.099 ; 
        end
        
        Turb.k_sigpsi = 2713;
        Turb.Window_size = 120;
    
    case 2 %% Axial impulse turbine    
        
    case 3
        %% Biradial turbine
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
        Turb.Om_vect = DB_reality.DB{1,1}.Vrota; %[100 110 110 110 120 130 140 150 165 170 190 195 215 220];
%         for i = 1:14
%             if Turb.Om_vect(i) > Turb.Om_co
%                 Turb.Om_vect(i) = Turb.Om_co*0.96;              
%             end
%         end 
        i = 0 ;
        while i < length(Turb.Om_vect)
            i = i+1 ;
            if Turb.Om_vect(i) > Turb.Om_co
                Turb.Om_vect(i) = Turb.Om_co*0.96;              
            end
        end 
        Turb.Om_opt = DB_reality.DB{1,1}.Vrota; %Ajouter le timing
            
        %n_CL = DB_reality.DB{1,1}.CL;
        % Variable speed
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

% Other constants           
WEC.p_star2psi = p_at / ( rho_at * Turb.D^2 );
WEC.c_at2 = gamma * p_at / rho_at;
WEC.GM1_G = ( gamma - 1.0 ) / gamma; 
end
