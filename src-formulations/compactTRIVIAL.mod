param a{M,N};

#param b := 10; #length of the hypercube containing the points
param b := 100; #length of the hypercube containing the points

param dUB := sqrt(n)*b;

var x{I}, binary;
var gamma{K};

var d{I} >= 0, <= dUB;

minimize objective:
  sum{(i,j) in I} d[i,j]^2;

subj to Assignment{i in M}:
  sum{j in K: (i,j) in I} x[i,j] = 1;

s.t. p_def1{(i,j) in I}:
  d[i,j] >= (sum{l in N} a[i,l]*w[j,l] - gamma[j]) - dUB*( 1 - x[i,j] );

s.t. p_def2{(i,j) in I}:
  d[i,j] >= -(sum{l in N} a[i,l]*w[j,l] - gamma[j]) - dUB*( 1 - x[i,j] );