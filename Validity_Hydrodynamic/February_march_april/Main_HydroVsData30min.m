clc
clearvars
close all

addpath Data

%%
%load DB_10min_april.mat
load DB_30min_feb_08h.mat
DB_trans = DB_30min_feb_08h; 

H0 = DB_trans.DB{1,1}.H0;
Pos_air = DB_trans.DB{1,1}.Position_air;
Rel_Pos_air = H0 - Pos_air;
Vitesse_column = DB_trans.DB{1,1}.Vitesse_column;
p_star_data = DB_trans.DB{1,1}.pPort/101325;
mean_H0 = mean(H0);

%%
% Simulation times
Times.Tsim = 1800; % Simulation duration
Times.Ts   = 0.05;   % Time step

% Initialise data from plant and turbine
% 1: Wells
% 2: Not used
% 3: Bi-radial
Turbine_type = 3;

% Choose control algo
% 0: Batch BeP
% 1: FS PI-ctrl
% 2: VS TorqueLaw k_opt
% 3: VS TorqueLaw best fit
% 4: MPC
CL = 1; 

if CL == 2
    a = 0.0011;
    b = 2;
    u_v = 1;
else 
    a =   0.0009241;%  (7.603e-05, 0.001772)
    b =       2.046;%
end

% Initialise waves
Time.ini = '2018-02-01 08:00:00.000';
Time.outi= '2018-02-01 08:30:00.000';
WaveData = Init_waves_Data_V2(Time);
Times.Trp  = 0.5 * WaveData.Tp;
[WEC, Turb] = Init_Mutriku_Data(DB_trans,Turbine_type,CL,mean_H0);

%%

tic
sim DATA_MTK_2.slx
toc

%% PLOT
t_data = linspace(Times.Ts,Times.Tsim,length(Pos_air));
t_calcule = linspace(Times.Ts,Times.Tsim,length(VelPos(:,1))); %t_calcule = Times.Ts : Times.Ts*2 : Times.Tsim ;

Rel_Pos_air = interp1(t_data,Rel_Pos_air,t_calcule);
Vitesse_column = interp1(t_data,Vitesse_column,t_calcule);
p_star_data = interp1(t_data,p_star_data,t_calcule);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position
size = 20;

figure()
subplot(2,1,1)
plot(t_calcule,Rel_Pos_air,t_calcule,VelPos(:,2))
legend('Mesure','Model')
title('Position')
set(gca,'fontsize',size)
subplot(2,1,2)
plot(t_calcule,Vitesse_column,t_calcule,VelPos(:,1))
legend('Mesure','Model')
title('Speed')
set(gca,'fontsize',size)

%%
figure()
% FFT1
Fs = 10;
L = length(Rel_Pos_air);
NFFT =  2^nextpow2(L);
outi = NFFT/20+1;
f = Fs/(2*pi)*linspace(0,1,outi);
fft_p_mdl = abs(fft(VelPos(:,2),NFFT)/L);
fft_p_meas  = abs(fft(Rel_Pos_air,NFFT)/L);

subplot(3,1,1)
plot(f,fft_p_meas(1:int32(outi)),f,fft_p_mdl(1:int32(outi)))
title('Spectrum of position')
legend('meas','model')
grid on
set(gca,'fontsize',size)

% FFT2
Fs = 10;
L = length(Vitesse_column);
NFFT =  2^nextpow2(L);
outi = NFFT/20+1;
f = Fs/(2*pi)*linspace(0,1,outi);
fft_p_mdl = abs(fft(VelPos(:,1),NFFT)/L);
fft_p_meas  = abs(fft(Vitesse_column,NFFT)/L);

subplot(3,1,2)
plot(f,fft_p_meas(1:int32(outi)),f,fft_p_mdl(1:int32(outi)))
title('Spectrum of speed')
legend('meas','model')
grid on
set(gca,'fontsize',size)

% FFT3
Fs = 10;
L = length(p_star_data);%length(Adim_res(:,4));
NFFT =  2^nextpow2(L);
outi = NFFT/20+1;
f = Fs/(2*pi)*linspace(0,1,outi);
fft_p_mdl = abs(fft(Adim_res(:,4),NFFT)/L);
fft_p_meas  = abs(fft(p_star_data,NFFT)/L);

subplot(3,1,3)
plot(f,fft_p_meas(1:int32(outi)),f,fft_p_mdl(1:int32(outi)))
legend('meas','model')
title('Spectrum of relavive pressure')
xlabel('frequence - rad/s')
grid on
set(gca,'fontsize',size)