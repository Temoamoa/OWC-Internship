%% DB access
clear 
close all
clc

Days = 1;
Month = '2-';
Hours = ["02:00:00" "02:30:00"];
%Time.ini = '2018-04-01 09:50:00.000';
%Time.outi= '2018-04-01 10:00:00.000';

DB_30min_feb_02h = DBaccess(Month,Days,Hours);
DB_30min_feb_02h = Mejorar(DB_30min_feb_02h,Days);
% save 'DataBase\DB_FullFev_without.mat' DB_FullFev_withoutstates
% save 'C:\Users\Temo\Desktop\Matlab_Mutriku_TM\Biradial\Data\DB_FullFev_without.mat' DB_FullFev_withoutstates
% save 'C:\Users\Temo\Desktop\Matlab_Mutriku_TM\Turbine Wells\Data\DB_FullFev_without.mat' DB_FullFev_withoutstates

save 'C:\Users\Temo\Desktop\Mutriku_W2W\Mutriku_NumModel_V5\Hydrodynamic\DB_30min_feb_02h.mat' DB_30min_feb_02h

function DB = DBaccess(Month,day,hrs)
NameFields = ["Time","state","CL","H0","Air level","VitesseColumn","Vitesse Rota Turbine ",...
        "Air pressure port site"];
    
Year =' "2018-'; 
Days = 1:31;

tic
for i = day
    fprintf('\n Now day %d',i);

    DB.Tlapse = strcat(Year,Month,num2str(Days(i))," ",char(hrs(1)),'" AND',Year,Month,num2str(Days(i))," ",char(hrs(2)),'"');
    %DB.Tlapse = strcat(' "2018-',num2str(i),'-26 09:30:00" AND "2018-',num2str(i),'-26 15:00:00" ');
    DB.Tlapse = char(DB.Tlapse);
    DB.Fields = ' Timestamp,Col017,Col002,Col023,Col011,Col026,Col052,Col010';

    DB.out{i,1} = getDBdataTM(DB.Tlapse,DB.Fields,NameFields);
end
fprintf('\n Results ready !\n');
toc

end

function structure = Mejorar(DB,days)

structure = struct();

for i = days %1 : size(DB.out,1)
    a = 1;
        
    DataJour = DB.out{i,a};

    if strcmp(DB.out{i,a},'No Data')
        structure.DB{i,a} = 'No Data';
    else
        Time = char(DataJour(2:end,1));
        State = str2num(char(DataJour(2:end,2)));
        CL = str2num(char(DataJour(2:end,3)));
%         M_T = str2num(char(DataJour(2:end,4)));
%         Torque = str2num(char(DataJour(2:end,3)));
        H0 = str2num(char(DataJour(2:end,4))); %Avg Air level
        IAS = str2num(char(DataJour(2:end,5))); %#ok<*ST2NM> %Air level PositionColonneAir Vc
        V_column = str2num(char(DataJour(2:end,6))); %Vitesse column of air
        %V_water = -V_volume;
        Vrota = str2num(char(DataJour(2:end,7))); %Rotation speed of the turbine
        pPort = str2num(char(DataJour(2:end,8)));

        %Relative Air and water Column
        %RIWS = H0 - IAS; %x dans la litérrature, quand la Interface Air Column descend la Interface Water Column monte 
        %RIAS = IAS - H0; %Relative IAS

        structure.DB{i,a}.Time = Time;
        structure.DB{i,a}.State = State;
        structure.DB{i,a}.CL = CL;
%         structure.DB{i,a}.M_T = M_T;
%         structure.DB{i,a}.Torque = Torque;
        structure.DB{i,a}.H0 = H0;
        structure.DB{i,a}.Position_air = IAS;
        %structure.DB{i,a}.Position_eau = RIWS;
        structure.DB{i,a}.Vitesse_column = V_column;
        structure.DB{i,a}.Vrota =  Vrota;
        structure.DB{i,a}.pPort =  pPort;
    end
end
end

%% DB Fieldname id
% Control Law Number        Col002	
% TestRunNumber             Col003	
% STA Current State         Col004	
% Vrel Sigma 60             Col005	m/s
% Vrel Sigma 300            Col006	m/s
% Pressure Sigma 60     	Col007	Pa
% Pressure Sigma 300        Col008	Pa
% Motor Torque              Col009	Nm
% Air pressure Port side	Col010          BONNNNNN
% Air Level               Col011	m Position water level Chamber
% Generator Temp 1          Col012	ºC
% Generator Temp 2          Col013	ºC
% Generator Temp 3          Col014	ºC
% Pressure RMS              Col015	Pa
% 	Col016	
% Current State             Col017	4 -> normal operation
% Run Time                  Col018	
% 	Col019	
% Damper Position           Col020	º
% Galeria Relative Humidity	Col021	%
% Sigma water velocity      Col022	m/s
% Avg Air Level           Col023	m H0
% Sigma water level         Col024	m
% Sigma water velocity      Col025	m/s
% Water velocity filt       Col026	m/s Vitesse !!!
% HSSV Position             Col027	mm
% HSSV Open                 Col028	1=Open/0=Closed
% Vibration Sensor          Col029	
% Flow                      Col030	m3/s
% 	Col031	
% Alarm Status              Col032	
% Pressure Sensor 1         Col033	Pa
% Pressure Sensor 2         Col034	Pa BONNNNNNNN
% Pressure Sensor 3         Col035	Pa
% Pressure Sensor 4         Col036	Pa
% Pressure Sensor 5         Col037	Pa 
% Pressure Sensor 6         Col038	Pa 
% Pressure Sensor 7         Col039	Pa
% Pressure Sensor 8         Col040	Pa
% Pressure Sensor 9         Col041	Pa
% Pressure Sensor 10        Col042	Pa BOOONNNNNNNNNNN
% Accelerometer 1           Col043	
% Accelerometer 2           Col044	
% Temperature Sensor K1     Col045	ºC
% Temperature Sensor K2     Col046	ºC
% Antifog                   Col047	
% PressureMean_300          Col048	Pa
% RT Wave elevation        	Col049	m
% RT Wave elevation no tide Col050	m
% Drive 1 Speed Ref         Col051	rad/s
% Drive 1 Speed Feedback	Col052	rad/s Vitesse rotation de la Turbine 
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
