clc
clear

load ('SS_RBRII.mat')
load ('Mutriku_SS.mat')
load ('Err_Allfev.mat')

% Year = SSRBRII(:,1);
% month = SSRBRII(:,2);
Day = SSRBRII(:,3);
Hour = SSRBRII(:,4);

Te = SSRBRII(:,6); % Période énergétique
Tp = SSRBRII(:,9); 
Hs = SSRBRII(:,12); % Hauteur

Indice_Fev = 12:683;%12:35 %+ (d-1)*24;
Indice_mars = 684:1427;
Indice_avril = 1427:1848;

Swell_fev = [Day(Indice_Fev) Hour(Indice_Fev) Hs(Indice_Fev) Te(Indice_Fev)]; 
Swell_mars = [Day(Indice_mars) Hour(Indice_mars) Hs(Indice_mars) Te(Indice_mars)];
Swell_avril = [Day(Indice_avril) Hour(Indice_avril) Hs(Indice_avril) Te(Indice_avril)]; %Tp(Indice_avril)

%%
Err_mean_poly = zeros(size(Err_Allfev_poly,1)/2,size(Err_Allfev_poly,2));
Err_mean_isen = zeros(size(Err_Allfev_poly,1)/2,size(Err_Allfev_poly,2));
Err_mean_poly_isen = zeros(size(Err_Allfev_poly,1)/2,size(Err_Allfev_poly,2));
for j = 1:size(Err_Allfev_poly,2)
    for i = 2:size(Err_Allfev_poly,1)
        if mod(i,2) == 0
            if Err_Allfev_poly(i-1,j) == 0 && Err_Allfev_poly(i,j) ~= 0
                Err_mean_poly(i/2,j) = Err_Allfev_poly(i,j);
                Err_mean_isen(i/2,j) = Err_Allfev_isen(i,j);
                Err_mean_poly_isen(i/2,j) = Err_isen_poly_Allfev(i,j);
            elseif Err_Allfev_poly(i-1,j) ~= 0 && Err_Allfev_poly(i,j) == 0
                Err_mean_poly(i/2,j) = Err_Allfev_poly(i-1,j);
                Err_mean_isen(i/2,j) = Err_Allfev_isen(i-1,j);
                Err_mean_poly_isen(i/2,j) = Err_isen_poly_Allfev(i-1,j);
            else
            Err_mean_poly(i/2,j) = (Err_Allfev_poly(i-1,j)+Err_Allfev_poly(i,j))/2;
            Err_mean_isen(i/2,j) = (Err_Allfev_isen(i-1,j)+Err_Allfev_isen(i,j))/2;
            Err_mean_poly_isen(i/2,j) = (Err_isen_poly_Allfev(i-1,j)+Err_isen_poly_Allfev(i,j))/2;
            end
        end
    end
end

