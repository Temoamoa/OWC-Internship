function DBout = getDBdataTM(Tlapse,Fields,NameFields)

    conn = database('OPERADB','TCN','Tecnalia2017');
    dbname = ' FROM operadb.t9data ';
    where  = strcat(' WHERE TimeStamp BETWEEN ', Tlapse);%,' AND (Col017=4) ');
    sqlquery = strcat('SELECT ', Fields, dbname, where,' ORDER BY id'); 
    curs = exec(conn,sqlquery);
    if curs.Message, curs.Message
    end
    curs = fetch(curs);
    
    if not(strcmp(curs.Data,'No Data'))
        a = [NameFields ; curs.Data];
        DBout =  a;
    else
        DBout =  curs.Data;
    end

end