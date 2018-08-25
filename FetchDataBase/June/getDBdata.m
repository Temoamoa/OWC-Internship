function DBout = getDBdata(Tlapse,Fields,CL)
    conn = database('OPERADB','TCN','Tecnalia2017');
    dbname = ' FROM operadb.t9data ';
    where  = strcat(' WHERE TimeStamp BETWEEN ', Tlapse, ' AND (Col017=4)AND (Col002=',num2str(CL),') ');
    sqlquery = strcat('SELECT ', Fields, dbname, where,' ORDER BY id'); 
    curs = exec(conn,sqlquery);
    if curs.Message, curs.Message
    end
    curs = fetch(curs);
    if not(strcmp(curs.Data,'No Data'))
    DBout =  curs.Data;
    end
end