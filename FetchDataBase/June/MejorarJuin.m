function structure = MejorarJuin(DB,Waves,Spect,Welev)

    structure = struct();
    DataJour = DB.out;

%     if strcmp(DB.out,'No Data')
%     structure.DB = 'No Data';
%     else
    Time = cell2mat(DataJour(:,1));
    State = cell2mat(DataJour(:,2));
    CL = cell2mat(DataJour(:,3));
    %         M_T = str2num(char(DataJour(2:end,4)));
    %         Torque = str2num(char(DataJour(2:end,3)));
    H0 = cell2mat(DataJour(:,4)); %Avg Air level
    IAS = cell2mat(DataJour(:,5)); %#ok<*ST2NM> %Air level PositionColonneAir Vc
    V_column = cell2mat(DataJour(:,6)); %Vitesse column of air
    %V_water = -V_volume;
    Vrota = cell2mat(DataJour(:,7)); %Rotation speed of the turbine
    pPort = cell2mat(DataJour(:,8));

    %Relative Air and water Column
    %RIWS = H0 - IAS; %x dans la litérrature, quand la Interface Air Column descend la Interface Water Column monte 
    %RIAS = IAS - H0; %Relative IAS
    
    i = 1;
    a = 1;
    
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
    structure.DB{i,a}.Waves = Waves;
    structure.DB{i,a}.Spect = Spect;
    structure.DB{i,a}.Welev = Welev;
    
end