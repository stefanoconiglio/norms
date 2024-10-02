#embedded antisymm: w[j,1] is forced to be >= 0

set N1 := N diff {1};
var w_p{K,N1 } >= 0, <= 1;
var w_n{K,N1 } >= 0, <= 1;
var s{K,N1} binary;

s.t. Norm1_1{j in K, l in N1}:
	w_p[j,l] - w_n[j,l] = w[j,l];
	
s.t. Norm1_2{j in K, l in N1}:
	w_p[j,l] <= s[j,l];	

s.t. Norm1_3{j in K, l in N1}:
	w_n[j,l] <= 1-s[j,l];	
	
s.t. Norm1_4{j in K}:
	w[j,1] + sum{l in N1} (w_p[j,l] + w_n[j,l]) >= 1;

