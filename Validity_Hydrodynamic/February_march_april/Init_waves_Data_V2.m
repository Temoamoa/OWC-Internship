%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Numercial modelling of Mutriku OWC                      %
%                                           OPERA Project                 %
%  F.X. Faÿ                                                               %
%  July 2016                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Waves = Init_waves_Data_V2(Time)
% Random generator starts always with the same seed
rng(17);
load Rad_SS.mat
% Vector of freqs
wstep =  Fint(3,1);           % Fundamental frequency discretization
wini =  Fint(3,1);
wend =  Fint(3,end);
wlen = wend - wini;
w_discr = Fint(3,:);% wini:wstep:wend;      % Frequency discretization

% Resource
% load Mutriku_SS.mat
% Hs = Mutriku(SS,1);
% Tp = Mutriku(SS,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hs, Tp from Ocean (Data)
load RBR_QC
Wave = getRessource_V2(Time,RBR_DB,Bimep);
Hs = Wave.SS(1);
Tp =  Wave.SS(2);
f = Wave.Spectre(:,1)'*2*pi;
Sxx = Wave.Spectre(:,2)';

w_nb = length(w_discr);

% Compute the perturbed intervals
dwr = zeros( 1, w_nb );
sum_dwr = 0.0;

for i = 1:w_nb
     if rand() >= 0.5
      sign = +0.2;
    else
      sign = -0.2;
    end

    dwr(i) = ( 1.0 + sign * rand() ) * wstep;
    sum_dwr = sum_dwr + dwr(i);
end

% scale the perturbed intervals to have the original wlen
% to get w_spec(end) == wend 
sum_dwr = sum_dwr - ( dwr(1) + dwr(end) ) / 2.0;
ratio = wlen * 0.999999 / sum_dwr; % avoid round up error

for i = 1:w_nb
    dwr(i) = dwr(i) * ratio;
end

%
w_discr = zeros( 1, w_nb ); 

for i = 1:w_nb
     if i == 1
        w_discr(1) = wini;           
    else
        w_discr(i) = w_discr(i-1) + 0.5 * ( dwr(i) + dwr(i-1) );
     end
end

% Generation of the wave spectrum 
%Sw2 = wave_spectrum_Data(Hs,Tp,w_discr);%*.35;

Sw = interp1(f, Sxx, w_discr);
i = 1;
while i < w_nb
    if isnan(Sw(i))
        Sw(i) = 0;
    end
    i = i + 1;
end

% figure()
% plot(w_discr/2.28,Sw,f,Sxx)
% legend('interp1','Sxx')

%Sw2 = wave_spectrum_Data(Hs,Tp,w_discr);%*.35;

% Sw3 = zeros(1,w_nb);
% Sw3(1:226) = Sxx;
% figure()
% plot(w_discr,Sw,w_discr,Sw3)
% legend('interp','non interp')

% Compute complexe Amplitudes and Phases
W_Am = zeros(1,length(w_discr));
W_Theta = zeros(1,length(w_discr));

for j = 1:length(w_discr)      
     W_Am(1,j) = (2*Sw(1,j)*dwr(j))^0.5; 
     W_Theta(1,j) = 2*pi*rand();
end
    
%% Output data
%Waves.SS = SS;
Waves.welev = Wave.welev;
Waves.Hs = Hs;
Waves.Tp = Tp;
Waves.Wam  = W_Am;
Waves.Wpha = W_Theta;
Waves.w_discr = w_discr;
Waves.w_nb = w_nb;
end


