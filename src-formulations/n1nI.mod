#embedded antisymm: in nI

var w_p{K,N} >= 0, <= 1;
var w_n{K,N} >= 0, <= 1;
var s{K,N} binary;

s.t. Norm1_1{j in K, l in N}:
	w_p[j,l] - w_n[j,l] = w[j,l];
	
s.t. Norm1_2{j in K, l in N}:
	w_p[j,l] <= s[j,l];	

s.t. Norm1_3{j in K, l in N}:
	w_n[j,l] <= 1-s[j,l];	
	
s.t. Norm1_4{j in K}:
        sum{l in N} (w_p[j,l] + w_n[j,l]) >= 1;


var u{K,N} binary;

s.t. Norminf_sqrtn_1{j in K, l in N}:
#     w[j,l] >= 1/sqrt(n)* ( 1 - (1+sqrt(n))*(1-u[j,l]) );
#     w[j,l] >= 1/sqrt(n) + (1-u[j,l]) * (-1 - 1/sqrt(n));
     w[j,l] >= -1 + u[j,l] * ( 1 + 1/sqrt(n) );
	
s.t. Norminf_sqrtn_2{j in K}:
	sum{l in N} u[j,l] = 1;
