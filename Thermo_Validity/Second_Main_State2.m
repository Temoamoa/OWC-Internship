clc
clearvars
close all


%load Data\DB_Full1Day02.mat
%load Data\DB_mars_SS4-5.mat
%load Data\DB_mars_SS9-10.mat
%load Data\DB_1Fev_without.mat
load DB_FullFev_without.mat
load Constante.mat
%load Data\DB_mars_5h
%%
DB_trans = DB_FullFev_withoutstates;%DB_mars;%DB_1Fev;%DB_Full1Day02;%DB_03_SS9_10;%DB_03_SS4_5;%DB_Full1Day02;%DB_Biradial02_1h;

Day = 1 : size(DB_trans.DB,1);
Demi_heure = 1 : size(DB_trans.DB,2); %47;

Err_Allfev_poly = zeros(size(DB_trans.DB,1),size(DB_trans.DB,2)); %28x48 pour l'instant
Err_Allfev_isen = zeros(size(DB_trans.DB,1),size(DB_trans.DB,2)); %28x48 pour l'instant

%% Check state of turbine
%%
% jour = size(DB_trans.DB,1);
% Demi_heure = size(DB_trans.DB,2);
% 
% j = 1;
% while j <= jour
%     
%     a = 1;
%     State_m = [];
%     
%     while a <= Demi_heure
%         if strcmp(DB_trans.DB{j, a},'No Data')
%             State_m = [State_m; zeros(7200,1)];
%         else
%             if ~isempty(DB_trans.DB{j, a})
%                 State_m = [State_m; DB_trans.DB{j, a}.State];
%             else
%                 State_m = [State_m; zeros(7200,1)];
%             end
%         end
%         a = a + 1;
%     end
%     
%     time = linspace(0,24,length(State_m));
%     plot(time,State_m)
%     title(sprintf('Feb Day %d',j))
%     ylabel('State')
%     xlabel('Time [h]')
%     drawnow
%     pause(1.5)
%     
%     j = j + 1;
% end
%%
path = "C:\Users\Temo\Desktop\Rapport\Stage MI4\Figures\OtherStyle\";
for j = 1%: size(DB_trans.DB,1) %Indice de jour
    for i = 1:5:21%size(DB_trans.DB,2)%29%26:36%1:18 %Demi_heure %15%1:20 %Apès y a des tests(plus state4) %Demi_heure %1:5 % Indice demi_heure
        if isempty(DB_trans.DB{Day(j),Demi_heure(i)})
            % Do nada
        elseif strcmp(DB_trans.DB{Day(j),Demi_heure(i)},'No Data')
            % No Data
        else 
            State = DB_trans.DB{Day(j),Demi_heure(i)}.State;
            V0 = DB_trans.DB{Day(j),Demi_heure(i)}.H0*Constante.S;
            Position_birad =DB_trans.DB{Day(j),Demi_heure(i)}.Position_air*Constante.S; %- DB_Biradial_trans.DB{Day,1}.H0*Constante.S;
            Vitesse_birad = DB_trans.DB{Day(j),Demi_heure(i)}.Vitesse_column*Constante.S;
            Omega = DB_trans.DB{Day(j),Demi_heure(i)}.Vrota;
            p_data = DB_trans.DB{Day(j),Demi_heure(i)}.pPort';
            p_star_data = p_data/Constante.p_atm;

            % Timing
            t0 = 0;
            tf_desire = 1800;

            % Obligatoirement entre 0 et 1800
            % Pour gérer les Indices temporels
            t = linspace(t0,tf_desire,length(DB_trans.DB{Day(j),Demi_heure(i)}.pPort))';
            t_ini = 1;
            t_end = length(t);
            Indice_temp = t_ini : t_end ;

            dt = t(2)-t(1); %incrément de temps
            dt_edo = dt;
            t_edo = t0:dt_edo:tf_desire;

            %% Résolution edo
            % CI
            Y0 =  p_star_data(Indice_temp(1)); 

            % Isen-Poly
            tic
            p_isen = f_gp_isen_State2(Y0,t,dt,State(Indice_temp),Omega(Indice_temp),Position_birad(Indice_temp),Vitesse_birad(Indice_temp),V0(Indice_temp))';
            p_poly = f_gp_poly_State2(Y0,t,dt,State(Indice_temp),Omega(Indice_temp),Position_birad(Indice_temp),Vitesse_birad(Indice_temp),V0(Indice_temp))';
            toc

            %% Err
            Err_Allfev_poly(j,i) = 100*(sum(abs(p_poly)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_poly)));
            Err_Allfev_isen(j,i) = 100*(sum(abs(p_isen)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_isen)));
            Err_halfHour = 100*(sum(abs(p_poly)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_poly)));
            %% Pour connaitre le déphasage des signaux
            [r, lag] = xcorr(p_poly,p_star_data(Indice_temp));
            [~,I] = max(abs(r));
            dcl = lag(I)*dt; %lag(I) %nbr d'élément à décalé
% 
%             %% Plot
%             figure()
%             subplot(3,1,1)
%             plot(t,State(Indice_temp))
%             title(sprintf('Err = %f',Err_1fev_poly(j,i)))
%             drawnow
% 
            size = 20;
            figure()
            subplot(2,1,1)
            plot(t_edo,p_poly,t+dcl,p_star_data(Indice_temp))
            title(sprintf('Model VS Data'))%, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
            xlabel('Time(s)')
            ylabel('p* [-]')
            legend('Model','Measure')
            set(gca,'fontsize',size)
            grid on
 
            % FFT
            Fs = 10;
            L = length(p_star_data(Indice_temp));
            NFFT =  2^nextpow2(L);
            outi = NFFT/20+1;
            f = Fs/(2*pi)*linspace(0,1,outi);
            fft_p_meas = abs(fft(p_star_data(Indice_temp),NFFT)/L);
            fft_p_mdl  = abs(fft(p_poly,NFFT)/L);

            subplot(2,1,2)
            plot(f,fft_p_mdl(1:int32(outi)-1),f,fft_p_meas(1:int32(outi)-1))
            %title('Spectrum of dimensionless relative pressure')
            title(sprintf('Spectrum of p*, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
            legend('Model','Measure')
            xlabel('frequence - rad/s')
            set(gca,'fontsize',size)
            %axis tight
            grid on
            drawnow
            name = strcat("Day_",num2str(Day(j)),"_half-hour_", num2str(Demi_heure(i)));
            Char = char(path+name);
            %saveas(gcf,Char,'png')
        end
    end
end

%%
Err_Allfev_poly = Err_Allfev_poly' ;
Err_Allfev_isen = Err_Allfev_isen' ;
Err_isen_poly_Allfev = Err_Allfev_isen - Err_Allfev_poly ;

for j = 1%: size(DB_trans.DB,1) %Indice de jour
    for i = 29%1:5:21%size(DB_trans.DB,2)%29%26:36%1:18 %Demi_heure %15%1:20 %Apès y a des tests(plus state4) %Demi_heure %1:5 % Indice demi_heure
        if isempty(DB_trans.DB{Day(j),Demi_heure(i)})
            % Do nada
        elseif strcmp(DB_trans.DB{Day(j),Demi_heure(i)},'No Data')
            % No Data
        else 
            State = DB_trans.DB{Day(j),Demi_heure(i)}.State;
            V0 = DB_trans.DB{Day(j),Demi_heure(i)}.H0*Constante.S;
            Position_birad =DB_trans.DB{Day(j),Demi_heure(i)}.Position_air*Constante.S; %- DB_Biradial_trans.DB{Day,1}.H0*Constante.S;
            Vitesse_birad = DB_trans.DB{Day(j),Demi_heure(i)}.Vitesse_column*Constante.S;
            Omega = DB_trans.DB{Day(j),Demi_heure(i)}.Vrota;
            p_data = DB_trans.DB{Day(j),Demi_heure(i)}.pPort';
            p_star_data = p_data/Constante.p_atm;

            % Timing
            t0 = 0;
            tf_desire = 1800;

            % Obligatoirement entre 0 et 1800
            % Pour gérer les Indices temporels
            t = linspace(t0,tf_desire,length(DB_trans.DB{Day(j),Demi_heure(i)}.pPort))';
            t_ini = 1;
            t_end = length(t);
            Indice_temp = t_ini : t_end ;

            dt = t(2)-t(1); %incrément de temps
            dt_edo = dt;
            t_edo = t0:dt_edo:tf_desire;

            %% Résolution edo
            % CI
            Y0 =  p_star_data(Indice_temp(1)); 

            % Isen-Poly
            tic
            p_isen = f_gp_isen_State2(Y0,t,dt,State(Indice_temp),Omega(Indice_temp),Position_birad(Indice_temp),Vitesse_birad(Indice_temp),V0(Indice_temp))';
            p_poly = f_gp_poly_State2(Y0,t,dt,State(Indice_temp),Omega(Indice_temp),Position_birad(Indice_temp),Vitesse_birad(Indice_temp),V0(Indice_temp))';
            toc

            %% Err
            Err_Allfev_poly(j,i) = 100*(sum(abs(p_poly)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_poly)));
            Err_Allfev_isen(j,i) = 100*(sum(abs(p_isen)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_isen)));
            Err_halfHour = 100*(sum(abs(p_poly)) - sum(abs(p_star_data(Indice_temp))))/(sum(abs(p_poly)));
            %% Pour connaitre le déphasage des signaux
            [r, lag] = xcorr(p_poly,p_star_data(Indice_temp));
            [~,I] = max(abs(r));
            dcl = lag(I)*dt; %lag(I) %nbr d'élément à décalé
% 
%             %% Plot
%             figure()
%             subplot(3,1,1)
%             plot(t,State(Indice_temp))
%             title(sprintf('Err = %f',Err_1fev_poly(j,i)))
%             drawnow
% 
            size = 20;
            figure()
            subplot(2,1,1)
            plot(t_edo,p_poly,t+dcl,p_star_data(Indice_temp))
            title(sprintf('Model VS Data'))%, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
            xlabel('Time(s)')
            ylabel('p* [-]')
            legend('Model','Measure')
            set(gca,'fontsize',size)
            grid on
 
            % FFT
            Fs = 10;
            L = length(p_star_data(Indice_temp));
            NFFT =  2^nextpow2(L);
            outi = NFFT/20+1;
            f = Fs/(2*pi)*linspace(0,1,outi);
            fft_p_meas = abs(fft(p_star_data(Indice_temp),NFFT)/L);
            fft_p_mdl  = abs(fft(p_poly,NFFT)/L);

            subplot(2,1,2)
            plot(f,fft_p_mdl(1:int32(outi)-1),f,fft_p_meas(1:int32(outi)-1))
            %title('Spectrum of dimensionless relative pressure')
            title(sprintf('Spectrum of p*, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
            legend('Model','Measure')
            xlabel('frequence - rad/s')
            set(gca,'fontsize',size)
            %axis tight
            grid on
            drawnow
            name = strcat("Day_",num2str(Day(j)),"_half-hour_", num2str(Demi_heure(i)));
            Char = char(path+name);
            %saveas(gcf,Char,'png')
        end
    end
end

%save Err_Allfev.mat Err_isen_poly_Allfev Err_Allfev_isen Err_Allfev_poly
%%
% figure()
% subplot(2,1,1)
% plot(t_edo,p_poly,t+dcl,p_star_data(Indice_temp))
% title(sprintf('Model VS Data, Err %2.1f%%',Err_halfHour))
% title('Offset rectified + Zoom')
% xlabel('Time(s)')
% ylabel('p* [-]')
% grid on
% 
% subplot(2,1,2)
% plot(t_edo,p_poly,t,p_star_data(Indice_temp))
% %title(sprintf('Model VS Data, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
% xlabel('Time(s)')
% title('Offset not rectified + Zoom')
% ylabel('p* [-]')
% grid on
% 
% % subplot(2,2,[3 4])
% % plot(f,fft_p_mdl(1:int32(outi)-1),f,fft_p_meas(1:int32(outi)-1))
% % %title('Spectrum of dimensionless relative pressure')
% % title(sprintf('Spectrum of p*'))%, Day %d half-hour %d, Err %2.1f%%',Day(j),Demi_heure(i),Err_halfHour))
% % legend('Model','Measure')
% % xlabel('frequence - rad/s')
% % %axis tight
% % grid on
% % drawnow
% % name = strcat("Day_",num2str(Day(j)),"_half-hour_", num2str(Demi_heure(i)),"Zoom2");
% % Char = char(path+name);
% % 
% % 
% % saveas(gcf,Char,'png')

%%
% Err_poly_hour = zeros(24,1);%length(Err_poly)/2);
% Err_isen_hour = zeros(24,1);
% for k = 1 : length(Err_mars_poly)
%     if mod(k,2) == 0
%     Err_poly_hour(k/2) = (Err_mars_poly(k)+Err_mars_poly(k+1))/2;
%     Err_isen_hour(k/2) = (Err_mars_isen(k)+Err_mars_isen(k+1))/2;
%     end
% end
% Err_isen_poly_hour = Err_isen_hour - Err_poly_hour;