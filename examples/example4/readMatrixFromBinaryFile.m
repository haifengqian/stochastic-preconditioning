function M = readMatrixFromBinaryFile( filename )

fid = fopen( filename );
dimension = fread( fid, 1, 'uint32' );
numnonzero = fread( fid, 1, 'uint32' );
rows = fread( fid, numnonzero, 'uint32' );
cols = fread( fid, numnonzero, 'uint32' );
values = fread( fid, numnonzero, 'double' );
fclose(fid);

M = sparse( rows+1, cols+1, values );
