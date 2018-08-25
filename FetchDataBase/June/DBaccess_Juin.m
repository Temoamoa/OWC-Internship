clear
clc

%% DB access
% clear all
% close all
% load CL1_TestTimeShort.mat
load AllCl_TestTime.mat
CL = 1;

%% 
tic
CL_Timevect = cell2mat(CL1_Timevect(:,1:2));
fprintf('\n Now CL %d',CL);

NameFields = ["Time","state","CL","H0","Air level","VitesseColumn","Vitesse Rota Turbine",...
        "Air pressure port site","WaveElev_ocean"];

for i=641 %632 %:521%length(CL_Timevect)
    fprintf('\n Now line: %d',i);
    Time.ini = CL_Timevect(i,1:19);
    Time.outi= CL_Timevect(i,20:38);

    DB_ini.Tlapse = strcat('"',Time.ini,'" AND "',Time.outi,'" ');
    DB_ini.Fields = ' Timestamp,Col017,Col002,Col023,Col011,Col026,Col052,Col034,Col049';  
    DB_ini.out = getDBdata(DB_ini.Tlapse,DB_ini.Fields,CL);
    
    % Mettre les noms au dessus c'est bête !!!
    %DB.out = getDBdata_NameFields(DB.Tlapse,DB.Fields,NameFields,CL); 
    
    [Waves,Spec,Welev_Rectifie] = getRessourceIsurki(Time,DB_ini.out(:,end));    
end

fprintf('\n Results ready !\n');
toc

Mejor_DB = MejorarJuin(DB_ini,Waves,Spec,Welev_Rectifie);
DB_12Juin_CL1_641 = Mejor_DB;

save 'C:\Users\Temo\Desktop\Mutriku_W2W\Mutriku_NumModel_V5\Hydrodynamic\DataBase\DB_12Juin_CL1_641.mat' DB_12Juin_CL1_641
save 'DB_12Juin_CL1_641.mat' DB_12Juin_CL1_641

% save CLxProdResults.mat ProdRes 
%% DB Fieldname id
% For D4.2 / Implementation 
% plots_DB(DB)

%%
% Control Law Number        Col002	
% TestRunNumber             Col003	
% STA Current State         Col004	
% Vrel Sigma 60             Col005	m/s
% Vrel Sigma 300            Col006	m/s
% Pressure Sigma 60     	Col007	Pa
% Pressure Sigma 300        Col008	Pa
% Motor Torque              Col009	Nm
% Air pressure Port side	Col010	
% Water Level               Col011	
% Generator Temp 1          Col012	ºC
% Generator Temp 2          Col013	ºC
% Generator Temp 3          Col014	ºC
% Pressure RMS              Col015	Pa
% 	Col016	
% Current State             Col017	
% Run Time                  Col018	
% 	Col019	
% Damper Position           Col020	º
% Galeria Relative Humidity	Col021	%
% Sigma water velocity      Col022	m/s
% Avg Water Level           Col023	m
% Sigma water level         Col024	m
% Sigma water velocity      Col025	m/s
% Water velocity filt	Col026	m/s
% HSSV Position             Col027	mm
% HSSV Open                 Col028	1=Open/0=Closed
% Vibration Sensor          Col029	
% Flow                      Col030	m3/s
% 	Col031	
% Alarm Status              Col032	
% Pressure Sensor 1         Col033	Pa
% Pressure Sensor 2         Col034	Pa
% Pressure Sensor 3         Col035	Pa
% Pressure Sensor 4         Col036	Pa
% Pressure Sensor 5         Col037	Pa
% Pressure Sensor 6         Col038	Pa
% Pressure Sensor 7         Col039	Pa
% Pressure Sensor 8         Col040	Pa
% Pressure Sensor 9         Col041	Pa
% Pressure Sensor 10        Col042	Pa
% Accelerometer 1           Col043	
% Accelerometer 2           Col044	
% Temperature Sensor K1     Col045	ºC
% Temperature Sensor K2     Col046	ºC
% Antifog                   Col047	
% PressureMean_300          Col048	Pa
% RT Wave elevation        	Col049	m
% RT Wave elevation no tide Col050	m
% Drive 1 Speed Ref         Col051	rad/s
% Drive 1 Speed Feedback	Col052	rad/s
% Drive 1 Total Current     Col053	A
% Drive 1 Active Current	Col054	A
% Drive 1 Torque Ref        Col055	%
% Drive 1 Out Hz            Col056	Hz
% Drive 1 Out V             Col057	V
% Drive 1 Out Power         Col058	W
% Drive 1 Out Power         Col059	W
% Drive 1 Out Power         Col060	W
% Regen Var Power           Col061	W
% Regen Total Current       Col062	A
% Regen Active Current      Col063	A
% Regen OutV                Col064	V
% Regen OutPower            Col065	W
% Regen BusDC               Col066	V
% Regen In Power1           Col067	VAR
% 	Col068	
% Speed Reference (IST)     Col069	Rpm
% Speed Feedback (IST)      Col070	Rpm
% Output Frequency (IST)	Col071	Hz
% Output Voltage (IST)      Col072	V
% Actual Current (IST)      Col073	A
% Actual Torque (IST)       Col074	N
% Active Power (IST)        Col075	W
