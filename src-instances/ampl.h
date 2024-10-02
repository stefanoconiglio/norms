#include"inclusions.h"

#ifndef __AMPL__
#define __AMPL__


/*!
 * \brief implements AMPL
 */
//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
class ampl : public algorithm
{
//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

//-------------
//! FIELDS
//-------------
private:
	int				K;				//! the maximum number of clusters

//-------------
//! CONSTRUCTOR
//-------------
public:
	//! constructor from problem_solution
	ampl(problem_solution * solution, int k = 1);

//-----------------------
//! COMPUTATIONAL SETTERS
//-----------------------
public:
	//! computes solution with Mangasarian's local search algorithm
	//! \return the total distance value
	pair<double,double> solve(const t_algorithm_choice ampl_choice, double ampl_bigM, double ampl_epsilon, bool use_continuous_relaxation);

private:
	//! write the .mod file for the monolithic model
	void generate_mod_monolithic_model(const t_algorithm_choice ampl_choice, bool use_continuous_relaxation = false);
	//! write the .dat file for the monolithic model
	void generate_dat_monolithic_model(const t_algorithm_choice ampl_choice, double ampl_bigM, double ampl_epsilon);
	//! write the .run file for the monolithic model
	void generate_run_monolithic_model();


//-------------
//! GETTERS
//-------------
public:


//--------------
//! PRINTERS
//--------------
public:


	void debug_on_video(double TS = 0, char* note = "nothing", bool complete_report = false, int level = 0);
	void debug_on_file(double TS = 0, char* note = "nothing", bool complete_report = false, int level = 0, file * debug_file = NULL);

};
#endif
