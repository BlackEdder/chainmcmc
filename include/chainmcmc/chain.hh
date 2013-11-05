/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of libmcmc.
  
 	libmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  libmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Gillespie. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef CHAIN_HH
#define CHAIN_HH

#include<vector>

#include <cppa/cppa.hpp>

namespace chainmcmc {
using namespace cppa;


/**
 * Parallel Metropolis Coupled Markov Chain Monte Carlo for Bayesian Phylogenetic Inference, 2004, Altekar et all.*/

/**
 * \brief Chain that runs mcmc process
 *
 * Chain is an actor and listens too certain commands (run no_generations
 * set temp, get temp)
 *
 * loglikelihood function, which takes parameters and data
 * parameters
 * priors, should always first do prior before calling loglikelihood, to see if it is usefull at all
 * temperature
 * data
 */

typedef double parameter_t;
typedef std::function<double( const std::vector<parameter_t> )> likelihood_t;
typedef std::function<double( const parameter_t )> prior_t;

namespace step {

	double mh_log_weight( const likelihood_t &ll,
			const std::vector<parameter_t> &pars,
			const std::vector<prior_t> &priors, const double &temperature );

	double mh_log_weight( const double &log_likelihood,
			const std::vector<parameter_t> &pars,
			const std::vector<prior_t> &priors, 
			const double &temperature );

	bool accept( std::mt19937 &eng, 
			const double &old_lweight,
			const double &new_lweight,
			const double &temperature = 1 );

	bool accept( std::mt19937 &eng, 
			const likelihood_t &ll, const double &old_ll,
			const std::vector<parameter_t> &old_pars,
			const std::vector<parameter_t> &new_pars,
			const std::vector<prior_t> &priors, double temperature = 1 );	

	bool accept( std::mt19937 &eng, const likelihood_t &ll, 
			const std::vector<parameter_t> &old_pars,
			const std::vector<parameter_t> &new_pars,
			const std::vector<prior_t> &priors, double temperature = 1 );

	parameter_t rparameter( std::mt19937 &eng, const parameter_t &par, 
			const double &sd );

	/**
	 * \brief Adapt the step size
	 *
	 * \param minimum_tries Only start adjusting after this number of tries
	 * \param min_level Adjust step size if acceptance rate below this level
	 * \param max_level Adjust step size if acceptance rate above this level
	 */
	double adapt_step_size( const double current_step_size,
			const size_t no_tries, const size_t no_accepts,
			const size_t minimum_tries, const double min_level, const double max_level );

	class ParameterState {
		public:
			size_t no_accepts = 0;
			size_t no_tries = 0;
			double sd = 0.001;
	};

	ParameterState adapt_parameter_sd( ParameterState && ps );

	/**
	 * \brief Stores the state of our mcmc chain
	 */
	class State {
		public:
			size_t generation = 0;
			double loglikelihood = NAN;
			std::vector<parameter_t> parameters;
			size_t current_parameter = 0; //! The parameter index we need to change
			std::vector<ParameterState> pss;
	};

	State step( std::mt19937 &eng, State && state, const likelihood_t &ll, 
			const std::vector<prior_t> &priors, 
			double temperature = 1 );
};

class Chain : public event_based_actor {
	public:
		Chain() {};
		/**
		 * \brief Start new chain
		 *
		 * \param loglikelihood function that calculates loglikelihood given a vector of parameters.
		 */
		Chain(  std::mt19937 &engine, 
		const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const std::vector<prior_t> &priors,
		const double & temperature = 1 );
		void init();

	protected:
		std::mt19937 rnd_engine;
		step::State state;
		double temperature = 1;
		const likelihood_t loglikelihood;
		std::vector<prior_t> priors;

		bool log_on = false;

		void step();
};

template<class T>
std::vector<T> fisherYatesKSubsets( std::vector<T> &v, 
		size_t k, std::mt19937 &rnd_engine ) { 
	size_t n = v.size()-1;
	for (size_t i=0; i<k; ++i) {
		std::uniform_int_distribution<size_t> distribution( 0, n );
		size_t j = distribution(rnd_engine);
		std::swap( v[j], v[n] );
		--n;
	}
	std::vector<T> ans(v.end()-k, v.end());
	return ans; 
};

/**
 * \brief Control the hot and cold chains
 *
 * Implemented following:
 * Parallel Metropolis Coupled Markov Chain Monte Carlo for Bayesian Phylogenetic Inference, 2004, Altekar et al
 *
 * General approach: 
 * Take a number of chains of different temperatures 
 * 	(one is the cold chain with temp 1).
 * Run them for a certain amount of steps
 * Choose two random ones and try to swap them
 * If accepted then swap temperatures.
 * 
 * For optimal performace we first choose the two to swap,
 * Then tell those two to run n generations, and all the others run 2n 
 * Then we try/perform the swap 
 * Tell the original two to run another n generations. Choose two and tell the
 * others to run to 3n generations etc
 */
class ChainController {
	public:
		ChainController( const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const std::vector<prior_t> &priors, size_t warm_up, size_t total_steps,
			size_t no_threads );

		ChainController( const likelihood_t &loglikelihood, 
		const std::vector<std::vector<parameter_t> > &pars_v,
		const std::vector<prior_t> &priors, size_t warm_up, size_t total_steps,
			size_t no_threads );


		void step();

	protected:
		std::mt19937 engine;
		size_t no_chains = 8;
		double dt = 0.1;
		std::vector<size_t> ids;
		std::map<size_t, actor_ptr> chains;

		size_t no_tries = 0;
		size_t no_accepts = 0;

		int warm_up = 0;

		size_t no_steps_between_swaps = 15; //! Try swap after this many steps
		void setup( const likelihood_t &loglikelihood, 
				const std::vector<std::vector<parameter_t> > &pars_v,
				const std::vector<prior_t> &priors,
				size_t no_chains );
		void run(const size_t total_steps);
};
};
#endif
