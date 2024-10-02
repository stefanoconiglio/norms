set M;
set N;
set K;
set I := setof{i in M, j in K: i >= j} (i,j);

param n := card(N);
