var u{K,N} binary;

s.t. Norminf_sqrtn_1{j in K, l in N}:
     w[j,l] >= -1 + u[j,l] * ( 1 + 1/sqrt(n) );
	
s.t. Norminf_sqrtn_2{j in K}:
	sum{l in N} u[j,l] = 1;
