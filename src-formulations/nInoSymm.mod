var u{K,N} binary;
var v{K,N} binary;

s.t. Norminf_sqrtn_1{j in K, l in N}:
     w[j,l] >= -1 + u[j,l] * ( 1 + 1/sqrt(n) );

s.t. Norminf_sqrtn_1_bis{j in K, l in N}:
     w[j,l] <= 1 + v[j,l] * ( -1 - 1/sqrt(n) );

s.t. Norminf_sqrtn_2{j in K}:
     sum{l in N} (u[j,l] + v[j,l]) = 1;
