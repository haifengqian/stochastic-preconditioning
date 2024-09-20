
system('./generate 50 m1_50_cube.matrix');
system('./generate 60 m2_60_cube.matrix');
system('./generate 70 m3_70_cube.matrix');
system('./generate 80 m4_80_cube.matrix');
system('./generate 90 m5_90_cube.matrix');
system('./generate 100 m6_100_cube.matrix');

system('./precondition m1_50_cube.matrix 0.35');
system('./precondition m2_60_cube.matrix 0.35');
system('./precondition m3_70_cube.matrix 0.35');
system('./precondition m4_80_cube.matrix 0.35');
system('./precondition m5_90_cube.matrix 0.35');
system('./precondition m6_100_cube.matrix 0.35');

global dimension
global matrix_nnz
global matrix_cond
global ICT_cholesky_nnz
global after_ICT_cond
global stoch_cholesky_nnz
global after_stoch_cond

runMatrix( 'm1_50_cube.matrix', 1 );
runMatrix( 'm2_60_cube.matrix', 2 );
plotCurve;
runMatrix( 'm3_70_cube.matrix', 3 );
close all;
plotCurve;
runMatrix( 'm4_80_cube.matrix', 4 );
close all;
plotCurve;
runMatrix( 'm5_90_cube.matrix', 5 );
close all;
plotCurve;
runMatrix( 'm6_100_cube.matrix', 6 );
close all;
plotCurve;
displayAllSummary;


