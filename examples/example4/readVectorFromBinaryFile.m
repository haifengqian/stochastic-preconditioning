function V = readVectorFromBinaryFile( filename )

fid = fopen( filename );
dimension = fread( fid, 1, 'uint32' );
V = fread( fid, dimension, 'uint32' );
V = V + 1;
fclose(fid);
