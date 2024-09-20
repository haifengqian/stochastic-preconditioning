function runMatrix( filename, index )

global dimension
global matrix_nnz
global matrix_cond
global ICT_cholesky_nnz
global after_ICT_cond
global stoch_cholesky_nnz
global after_stoch_cond

clear Q R
global Q
global R

Q = readMatrixFromTextFile( filename );
dimension( index ) = size(Q,1);
matrix_nnz( index ) = nnz( Q );
opts.issym = 1;
opts.maxit = 1000;
large_Q = eigs( Q, 1, 'la', opts);
small_Q = eigs( Q, 1, 'sa', opts);
matrix_cond( index ) = large_Q/small_Q;

p = symamd(Q);
Q = Q(p,p);
R = cholinc(Q,0.004);
ICT_cholesky_nnz( index ) = nnz( R );
large_precondQ = eigs( @precondFunc, size(Q,1), 1, 'la', opts );
small_precondQ = eigs( @precondFunc, size(Q,1), 1, 'sa', opts );
after_ICT_cond( index ) = large_precondQ/small_precondQ;

clear Q R p
global Q
global R

Q = readMatrixFromTextFile( filename );
p = readVectorFromBinaryFile( cat( 2, filename, '.ordering' ) );
Q = Q(p,p);
R = readMatrixFromBinaryFile( cat( 2, filename, '.cholesky' ) );
stoch_cholesky_nnz( index ) = nnz( R );
large_precondQ = eigs( @precondFunc, size(Q,1), 1, 'la', opts );
small_precondQ = eigs( @precondFunc, size(Q,1), 1, 'sa', opts );
after_stoch_cond( index ) = large_precondQ/small_precondQ;

clear Q R

displaySummary( index );

