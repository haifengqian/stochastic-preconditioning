function displaySummary( index )

global dimension
global matrix_nnz
global matrix_cond
global ICT_cholesky_nnz
global after_ICT_cond
global stoch_cholesky_nnz
global after_stoch_cond

fprintf('----------------- the %d-th matrix -----------------\n', index);
fprintf('dimension of matrix: %d\n', dimension(index));
fprintf('number of non-zero of the original matrix (diagonal + upper + lower): %d\n', matrix_nnz(index));
fprintf('condition number of the original matrix: %e\n', matrix_cond(index));
fprintf('number of non-zero of the ICT preconditioner (diagonal + upper): %d\n', ICT_cholesky_nnz(index));
fprintf('condition number after ICT preconditioning: %e\n', after_ICT_cond(index));
fprintf('number of non-zero of the stochastic preconditioner (diagonal + upper): %d\n', stoch_cholesky_nnz(index));
fprintf('condition number after stochastic preconditioning: %e\n', after_stoch_cond(index));
fprintf('---------------------------------------------------\n');

