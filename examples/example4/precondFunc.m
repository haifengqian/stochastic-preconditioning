function y = precondFunc(x)
global Q
global R

y = ( R' ) \ ( Q * ( R \ x ) );
