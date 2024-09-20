function M = readMatrixFromTextFile( filename )

fid = fopen( filename );
tline = fgets(fid);
dimension = sscanf(tline,'%d');
tline = fgets(fid);
entrynumber = sscanf(tline,'%d');
readbuffer = fscanf(fid,'%d %d %f',3*entrynumber);
fclose(fid);

M = sparse(readbuffer(1:3:3*entrynumber-2)+1,readbuffer(2:3:3*entrynumber-1)+1,readbuffer(3:3:3*entrynumber));
clear readbuffer
