function structure = MejorarOut(DB,days,hrs)

structure = struct();

for i = days %1 : size(DB.out,1)
    for a = hrs %1 : size(DB.out,2)
        fprintf('\n Hour %f',a/2);
        DataJour = DB.out{i,a};

        if strcmp(DB.out{i,a},'No Data')
            structure.DB{i,a} = 'No Data';
        else
            Time = char(DataJour(2:end,1));
            State = str2num(char(DataJour(2:end,2)));
            H0 = str2num(char(DataJour(2:end,3))); %Avg Air level
            IAS = str2num(char(DataJour(2:end,4))); %#ok<*ST2NM> %Air level PositionColonneAir Vc
            V_column = str2num(char(DataJour(2:end,5))); %Vitesse column of air
            %V_water = -V_volume;
            Vrota = str2num(char(DataJour(2:end,6))); %Rotation speed of the turbine
            pPort = str2num(char(DataJour(2:end,7)));

            %Relative Air and water Column
            %RIWS = H0 - IAS; %x dans la litérrature, quand la Interface Air Column descend la Interface Water Column monte 
            %RIAS = IAS - H0; %Relative IAS

            structure.DB{i,a}.Time = Time;
            structure.DB{i,a}.State = State;
            structure.DB{i,a}.H0 = H0;
            structure.DB{i,a}.Position_air = IAS;
            %structure.DB{i,a}.Position_eau = RIWS;
            structure.DB{i,a}.Vitesse_column = V_column;
            structure.DB{i,a}.Vrota =  Vrota;
            structure.DB{i,a}.pPort =  pPort;
        end
    end
end

end