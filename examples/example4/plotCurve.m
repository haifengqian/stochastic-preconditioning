function plotCurve()

global dimension
global after_ICT_cond
global after_stoch_cond

plot( dimension, after_ICT_cond, 'r' );
hold on
plot( dimension, after_stoch_cond, 'b' );
xlabel( 'Matrix Dimension' );
ylabel( 'Condition Number after Preconditioning' );
title( 'Condition Number Comparison between Equal-sized ICT and Stochastic Preconditioning. Red: ICT. Blue: Stochastic.' );

