#include"ampl.h"

void
ampl::solve_2norm_MINLP(Problem_instance & pi)
{
	ofstream ampl_model("MINLP.mod");
	
	ampl_model << "\
	set M;	#points							\n\
	set N;	#space dimensions						\n\
	set K;	#cluster space						\n\
	var x {M, K} binary;							\n\
	var gamma {K};								\n\
	var w {K, N} >= -1, <= 1;						\n\
	param a {M, N};								\n\
	\n\
	minimize obj:								\n\
	\n\
	sum{i in M, j in K} sum{l in N} (a[i,l]*w[j,l] - gamma[j])^2 * x[i,j];	\n\
	subject to									\n\
	\n\
	covering{i in M}:								\n\
	sum{j in K} ( x[i,j] ) = 1;						\n\
	\n\
	norm_2{j in K}:								\n\
	sum{l in N} w[j,l]^2 = 1;						\n\
	" << endl;

	ofstream ampl_data("MINLP.dat");
	pi.output_standard_ampl_data(ampl_data);

	ofstream ampl_runner("MINLP.run");

	#ifndef __ORLAB__
		ampl_runner << "option solver cplex;" << endl;
	#else
		ampl_runner << "option solver cplexamp100.3;" << endl;
	#endif
	ampl_runner << "\
	option cplex_options 'timing 1';							\n\
	option cplex_options $cplex_options 'integrality = 1e-20';			\n\
	model ampl_monolithic.mod;								\n\
	data ampl_monolithic.dat;								\n\
	solve;											\n\
	display obj > MINLP.output;								\n\
	display x > MINLP.output;								\n\
	display w > MINLP.output;								\n\
	display gamma > MINLP.output;								\n\
	" << endl;

	system("ampl MINLP.run");
}


void
ampl::solve_inf_norm_MIQP(Problem_instance & pi)
{
	ofstream ampl_model("inf_MIQP.mod");
	
	ampl_model << "\
	set M;	#points									\n\
	set N;	#space dimensions								\n\
	set K;	#cluster space								\n\
	var x {M, K} binary									\n\
	var d {M, K} >= 0;									\n\
	var gamma {K};										\n\
	var w {K, N} >= -1, <= 1;								\n\
	var tplus {K, N} >= -1, <= 1;								\n\
	var tminus {K, N} >= -1, <= 1;							\n\
	var z {K, N} binary;									\n\
	param a {M, N};										\n\
	param M = 1;										\n\
													\n\
	minimize obj:										\n\
													\n\
	sum{i in M} d[i]^2 									\n\
	subject to											\n\
													\n\
	covering{i in M}:										\n\
	sum{j in K} ( x[i,j] ) = 1;								\n\
													\n\
	dist_1 {i in M, j in K}:								\n\
	d[i] >= sum{l in D} ( a[i,l] * w[j,l] ) - gamma[j] - M * (1-x[i,j]);	\n\
													\n\
	dist_2 {i in M, j in K}:								\n\
	d[i] >= sum{l in D} ( - a[i,l] * w[j,l] ) + gamma[j] - M * (1-x[i,j]);	\n\
													\n\
	1_norm {j in K}:										\n\
	sum{l in N} tplus[j,l] + tminus[j,l] = 1; 					\n\
													\n\
	1_norm_tplus_minus {j in K, l in N}:						\n\
	tplus[j,l] - tminus[j,l] = w[j,l];							\n\
													\n\
	1_norm_activation{j in K, l in N}:							\n\
	tplus[j,l] <= z[j,l];									\n\
	tminus[j,l] <= 1-z[j,l];								\n\
	" << endl;

	ofstream ampl_data("inf_MIQP.dat");
	pi.output_standard_ampl_data(ampl_data);

	ofstream ampl_runner("inf_MIQP.run");

	#ifndef __ORLAB__
		ampl_runner << "option solver cplex;" << endl;
	#else
		ampl_runner << "option solver cplexamp100.3;" << endl;
	#endif
	ampl_runner << "\
	option cplex_options 'timing 1';							\n\
	option cplex_options $cplex_options 'integrality = 1e-20';			\n\
	model ampl_monolithic.mod;								\n\
	data ampl_monolithic.dat;								\n\
	solve;											\n\
	display obj > inf_MIQP.output;							\n\
	display x > inf_MIQP.output;								\n\
	display d > inf_MIQP.output;								\n\
	display w > inf_MIQP.output;								\n\
	display gamma > inf_MIQP.output;							\n\
	display z > inf_MIQP.output;								\n\
	display tplus > inf_MIQP.output;							\n\
	display tminus > inf_MIQP.output;							\n\
	" << endl;

	system("ampl inf_MIQP.run");
}

void
ampl::solve_1norm_MIQP(Problem_instance & pi)
{
	ofstream ampl_model("1MIQP.mod");
	
	ampl_model << "\
	set M;	#points									\n\
	set N;	#space dimensions								\n\
	set K;	#cluster space								\n\
	var x {M, K} binary									\n\
	var d {M, K} >= 0;									\n\
	var gamma {K};										\n\
	var w {K, N} >= -1, <= 1;								\n\
	var z {K, N} binary;									\n\
	param a {M, N};										\n\
	param M = 1;										\n\
	\n\
	minimize obj:										\n\
	\n\
	sum{i in M} d[i]^2 									\n\
	subject to											\n\
	\n\
	covering{i in M}:										\n\
	sum{j in K} ( x[i,j] ) = 1;								\n\
	\n\
	dist_1 {i in M, j in K}:								\n\
	d[i] >= sum{l in D} ( a[i,l] * w[j,l] ) - gamma[j] - M * (1-x[i,j]);	\n\
	\n\
	dist_2 {i in M, j in K}:								\n\
	d[i] >= sum{l in D} ( - a[i,l] * w[j,l] ) + gamma[j] - M * (1-x[i,j]);	\n\
	\n\
	infinite_norm {j in K, l in N}:							\n\
	w[j,l] >= 1 - 2*(1-z[j,l]);								\n\
	\n\
	infinite_norm_2{j in K}:								\n\
	sum{l in N} (z[j,l]) = 1;								\n\
	" << endl;

	ofstream ampl_data("1MIQP.dat");
	pi.output_standard_ampl_data(ampl_data);

	ofstream ampl_runner("1MIQP.run");

	#ifndef __ORLAB__
		ampl_runner << "option solver cplex;" << endl;
	#else
		ampl_runner << "option solver cplexamp100.3;" << endl;
	#endif
	ampl_runner << "\
	option cplex_options 'timing 1';							\n\
	option cplex_options $cplex_options 'integrality = 1e-20';			\n\
	model ampl_monolithic.mod;								\n\
	data ampl_monolithic.dat;								\n\
	solve;											\n\
	display obj > 1MIQP.output;								\n\
	display x > 1MIQP.output;								\n\
	display d > 1MIQP.output;								\n\
	display w > 1MIQP.output;								\n\
	display gamma > 1MIQP.output;								\n\
	display z > 1MIQP.output;								\n\
	" << endl;

	system("ampl 1MIQP.run");
}


//! reminder on how to fix variables
// 	ampl_model.write_line("fix x[1,2] := 1;");

