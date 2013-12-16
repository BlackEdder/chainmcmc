/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of chainmcmc.
  
  Chainmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Chainmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with chainmcmc. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#include "chainmcmc/chain.hh"

namespace chainmcmc {

/**
 * \brief step contains many small functions that are used for doing steps
 *
 * Design decisions: These functions are all self contained to improve testability. This design is influenced by functional programming paradigms, again mainly to improve testability and reusability.
 */
namespace step {
	double mh_log_weight( const likelihood_t &ll,
			const std::vector<parameter_t> &pars,
			const joint_prior_t &joint_prior, const double &temperature ) {
		double pr = joint_prior( pars );
		if ( pr == 0 )
			return log(0);
		return (log(pr) + ll( pars ))/temperature;
	}

	double mh_log_weight( const double &log_likelihood,
			const std::vector<parameter_t> &pars,
			const joint_prior_t &joint_prior, 
			const double &temperature ) {
		if (!std::isfinite( log_likelihood ))
			return log_likelihood;
		double pr = joint_prior( pars );
		if ( pr == 0 )
			return log(0);
		return (log_likelihood+log(pr))/temperature;
	}

	/**
	 * \brief Metropolis Hastings acceptance rule.
	 *
	 * \param temperature is actually unused here and assumed to be already included in the lweight parameters
	 */
	bool accept( std::mt19937 &eng, 
			const double &old_lweight,
			const double &new_lweight,
			const double &temperature ) {
		if (std::isnan( old_lweight))
			return true;
		if (std::isnan( new_lweight))
			return false;
		if (!std::isfinite( old_lweight ) && old_lweight > 0) 
			std::cerr << "Infinite likelihood detected. Probably numerical error." 
				<< std::endl;
		double log_ratio = new_lweight - old_lweight;
		if (log_ratio >= 0)
			return true;
		else {
			std::uniform_real_distribution<double> unif(0,1);
			return log(unif(eng)) < log_ratio;
		}
		return false;
	}

	bool accept( std::mt19937 &eng, 
			const likelihood_t &ll, const double &old_ll,
			const std::vector<parameter_t> &old_pars,
			const std::vector<parameter_t> &new_pars,
			const joint_prior_t &joint_prior, double temperature ) {

		double old_lweight = mh_log_weight( old_ll, old_pars, joint_prior, temperature );
		if (std::isnan( old_lweight ))
			return true;
		return accept( eng, old_lweight,
				mh_log_weight( ll, new_pars, joint_prior, temperature ), temperature );
	}

	bool accept( std::mt19937 &eng, const likelihood_t &ll, 
			const std::vector<parameter_t> &old_pars,
			const std::vector<parameter_t> &new_pars,
			const joint_prior_t &joint_prior, double temperature ) {

		// First check if none of the priors are zero before we call any likelyhood
		// function
		double pr = joint_prior( new_pars );
		if ( pr == 0 )
			return false;
		return accept( eng, ll, ll( old_pars ), old_pars, new_pars, 
				joint_prior, temperature );	
	}

	parameter_t rparameter( std::mt19937 &eng, const parameter_t &par, 
			const double &sd ) {
		std::normal_distribution<double> rnorm( 0, sd );
		if ( std::isfinite( par ) )
			return par + rnorm( eng );
		else if (par<0) // -inf
			return std::numeric_limits<double>::min()*100 + rnorm( eng );
		else
			return std::numeric_limits<double>::max()/100.0 + rnorm( eng );
	}

	double adapt_step_size( const double current_step_size,
			const size_t no_tries, const size_t no_accepts,
			const size_t minimum_tries, const double min_level, const double max_level,
			const size_t interval ) {
		double step_size = current_step_size;
		if (no_tries>minimum_tries && no_tries%interval == 0 ) {
			double alpha = ((double) no_accepts)/no_tries;
			if ( alpha > max_level)
				step_size*=1.05;
			else if ( alpha < min_level)
				step_size/=1.05;
		}
		if (step_size <= 0)
			step_size = std::numeric_limits<double>::min()*1000;
		return step_size;
	}


	ParameterState adapt_parameter_sd( ParameterState && ps ) {
		ps.sd = adapt_step_size( ps.sd, ps.no_tries, ps.no_accepts,
				100, 0.20, 0.25, 10 );
		return ps;
	}

	State step( std::mt19937 &eng, State && state, const likelihood_t &ll, 
			const joint_prior_t &joint_prior, bool adapting,
			double temperature ) {
		++state.generation;
		++state.current_parameter;
		state.current_parameter = state.current_parameter%state.parameters.size();
	
		bool accepted = true;
		auto proposed_parameters = state.parameters;

		proposed_parameters[state.current_parameter] = 
			rparameter( eng, proposed_parameters[state.current_parameter], 
					state.pss[state.current_parameter].sd	);
		++state.pss[state.current_parameter].no_tries;

		if (joint_prior( proposed_parameters ) == 0) {
			accepted = false;
		}

		double old_lweight;
		double proposed_ll;

		if (accepted) {
			if (std::isnan( state.loglikelihood )) {
				proposed_ll = ll( proposed_parameters );
			} else {
				proposed_ll = ll( proposed_parameters );
				old_lweight = mh_log_weight( state.loglikelihood, state.parameters,
						joint_prior, temperature );
				if (!std::isnan( old_lweight)) {
					double proposed_lweight = mh_log_weight( proposed_ll, proposed_parameters,
							joint_prior, temperature );
					accepted = accept( eng, old_lweight, proposed_lweight, temperature );
				}
			}
		}

		if (accepted) {
			++state.pss[state.current_parameter].no_accepts;

			state.loglikelihood = proposed_ll; 
			state.parameters = proposed_parameters;
		}
		if (adapting) {
			state.pss[state.current_parameter] =
				adapt_parameter_sd( std::move( state.pss[state.current_parameter] ) );
		}

		return state;
	}
};

Chain::Chain( std::mt19937 &engine, 
		const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const std::vector<prior_t> &priors,
		const double & temperature )
: rnd_engine( engine ),
	temperature( temperature ),
	loglikelihood( loglikelihood ),
	joint_prior( priors )
{
	state.parameters = parameters;
	for ( size_t i = 0; i < parameters.size(); ++i )
		state.pss.push_back( step::ParameterState() );
	state.pss[0].sd = 100;
	state.pss[1].sd = 1;
	state.pss[2].sd = 30000;
}

Chain::Chain( std::mt19937 &engine, 
		const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const joint_prior_t &joint_prior,
		const double & temperature )
: rnd_engine( engine ),
	temperature( temperature ),
	loglikelihood( loglikelihood ),
	joint_prior( joint_prior )
{
	state.parameters = parameters;
	for ( size_t i = 0; i < parameters.size(); ++i )
		state.pss.push_back( step::ParameterState() );
	state.pss[0].sd = 100;
	state.pss[1].sd = 1;
	state.pss[2].sd = 30000;
}
void Chain::init()  {
	become(
		on( atom( "run" ), arg_match ) >> [this]( const int no ) {
			for (size_t i = 0; i < (size_t)no; ++i) {
				step();
			}
		},
		on( atom( "run" ), arg_match ) >> [this]( 
			const int no, const bool start_logging ) {
			log_on = start_logging;
			for (size_t i = 0; i < (size_t)no; ++i) {
				step();
			}
		},
		on( atom( "run" ), arg_match ) >> [this]( const size_t no ) {
			for (size_t i = 0; i < no; ++i) {
				step();
			}
		},
		on( atom( "run" ), arg_match ) >> [this]( 
			const size_t no, const bool start_logging ) {
			log_on = start_logging;
			for (size_t i = 0; i < no; ++i) {
				step();
			}
		},
		on( atom( "no_adapt" ) ) >> [this]() {
			/*std::cout << "Sd ";
			for ( auto & ps : state.pss )
				std::cout << " " << ps.sd;
			std::cout << std::endl;*/
			adapting = false;
		},
		on( atom( "step" ) ) >> [this]() {
			step();
		},
		on( atom("log_weight") ) >> [this]() {
			// Returns the cold weight
			return step::mh_log_weight( state.loglikelihood, state.parameters,
				joint_prior, 1 );
		},
		on( atom("temp" ), arg_match ) >> [this]( const double &new_temp ) {
			temperature = new_temp;
		},
		on( atom("temp" ) ) >> [this]() {
			return temperature;
		},
		on( atom( "logger" ), arg_match ) >> [this]( 
				const actor_ptr &new_logger ) {
			usleep( 10000 );
			send( logger, atom("close") );
			usleep( 10000 );
			logger = new_logger;
		},
		on( atom("close" ) ) >> [this]() {
			self->quit();
			return atom("closed");
		}
	);
}

void Chain::step()  {
	state = step::step( rnd_engine, std::move( state ), loglikelihood, joint_prior,
			adapting, temperature );
	if (state.generation%50 == 0 && log_on) {
		std::stringstream s; // Collect output in stringstream for thread safety
		bool first = true;
		/*std::cout << state.loglikelihood << ", " << temperature << 
			", " << state.current_parameter << ": ";*/
		for ( auto & par : state.parameters ) {
			if (first) {
				s << par;
				first = false;
			} else
				s	<< "\t" << par;
		}
		send( logger, atom("append"), s.str() ); 
	}
}

ChainController::ChainController( const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const joint_prior_t &joint_prior, size_t warm_up, size_t total_steps,
		size_t no_chains, std::ostream &out ) 
	: no_chains( no_chains ), warm_up( warm_up ) {
		setup( loglikelihood, {parameters}, joint_prior, no_chains, out );

		run( total_steps );
	}


ChainController::ChainController( const likelihood_t &loglikelihood, 
		const std::vector<std::vector<parameter_t> > &pars_v,
		const joint_prior_t &joint_prior, size_t warm_up, size_t total_steps,
		size_t no_chains, std::ostream &out ) : no_chains( no_chains ), warm_up( warm_up ) {

	setup( loglikelihood, pars_v, joint_prior, no_chains, out );

	run( total_steps );
}

void ChainController::step() {
	warm_up -= no_steps_between_swaps;
	if (chains.size()>1) {
		fisherYatesKSubsets( ids, 2, engine );
		for ( size_t i = 2; i < ids.size(); ++i ) {
			bool log_on = false;
			if (warm_up<0 && ids[i]==0)
				log_on = true;
			send( chains[ids[i]], atom("run"), no_steps_between_swaps, log_on );
		}

		double weight1, weight2;
		double temp1, temp2;

		send( chains[ids[0]], atom("log_weight") );
		receive( 
			on(arg_match) >> [&weight1]( const double &ans ) {
				weight1 = ans;
			}
		);
		send( chains[ids[1]], atom("log_weight") );
		receive( 
			on(arg_match) >> [&weight2]( const double &ans ) {
				weight2 = ans;
			}
		);

		send( chains[ids[0]], atom("temp") );
		receive( 
			on(arg_match) >> [&temp1]( const double &ans ) {
				temp1 = ans;
			}
		);
		send( chains[ids[1]], atom("temp") );
		receive( 
			on(arg_match) >> [&temp2]( const double &ans ) {
				temp2 = ans;
			}
		);
		// Try swap
		// ratio = weight2^(1/temp1)*weight1^(1/temp2)/(weight2^(1/temp2)*weight1^(1/temp1))
		double log_ratio = ((temp2-temp1)*weight2+(temp1-temp2)*weight1)/(temp1*temp2);
		bool accept = false;
		if (std::isfinite(log_ratio)) { // Only consider switching when ratio is a finite number
			if (log_ratio > 0)
				accept = true;
			else {
				std::uniform_real_distribution<double> unif(0,1);
				if (log(unif(engine))<log_ratio)
					accept = true;
			}
		}
		++no_tries;
		if (accept) {
			++no_accepts;

			dt = 1.0/step::adapt_step_size( 1.0/dt, no_tries, no_accepts,
				10, 0.2, 0.6 );

			send( chains[ids[0]], atom("temp"), (1+dt*ids[1]) );
			send( chains[ids[1]], atom("temp"), (1+dt*ids[0]) );
		}

		bool log_on = false;
		if (warm_up<0 && ids[1]==0)
			log_on = true;

		send( chains[ids[1]], atom("run"), no_steps_between_swaps, log_on );
	}

	bool log_on = false;
	if (warm_up<0 && ids[0]==0)
		log_on = true;
	send( chains[ids[0]], atom("run"), no_steps_between_swaps, log_on );
}

	void ChainController::setup( const likelihood_t &loglikelihood, 
			const std::vector<std::vector<parameter_t> > &pars_v,
			const joint_prior_t &joint_prior,
			size_t no_chains, std::ostream &out ) {
		for ( size_t i = 0; i < no_chains; ++i ) {
			std::mt19937 eng;
			eng.seed( engine() );
			ids.push_back( i );
			double temp = (1+dt*i);
			auto chain = spawn<Chain>( eng, loglikelihood, pars_v[i%pars_v.size()],
					joint_prior, temp );
			logger = spawn<Logger>( out );
			send( chain, atom("logger"), logger );
			send( chain, atom("run"), no_steps_between_swaps );
			chains[i] = chain;
		} 
	}
	
void ChainController::run( const size_t total_steps ) {
	size_t steps = (warm_up+total_steps)/no_steps_between_swaps;
	for ( size_t i = 0; i < steps; ++i )
		step();
	for ( auto & id_chain : chains ) 
		send( id_chain.second, atom("close") );
	size_t i = 0; 
	receive_for( i, chains.size() ) (
		on( atom("closed") ) >> []() {}
	);
	send( logger, atom("close") );
	receive(
		on( atom("closed") ) >> []() {}
	);
}

FPChainController::FPChainController( const likelihood_t &loglikelihood, 
		const std::vector<std::vector<parameter_t> > &pars_v,
		const joint_prior_t &joint_prior, size_t warm_up, size_t total_steps,
		std::ostream &out ) : 
			warm_up( warm_up ), total_steps( total_steps ), out( out ) {
		setup( loglikelihood, pars_v, joint_prior );
	}

	likelihood_t FPChainController::heated_loglikelihood( 
			const double &temp, const likelihood_t &loglikelihood ) {
		return [temp, loglikelihood]( const trace::sample_t &pars ) {
			if (temp == 0)
				return 0.0;
			return temp*loglikelihood( pars );
		};
	}

	void FPChainController::setup( const likelihood_t &loglikelihood, 
			const std::vector<std::vector<parameter_t> > &pars_v,
			const joint_prior_t &joint_prior ) {
		for (size_t i = 0; i < n; ++i) {
			double temp = pow( ((double) i)/(n-1), c );
			traces[temp] = std::vector<trace::sample_t>();
			trace_actors[temp] = spawn<TraceLogger>( traces[temp] );
			FPChainState state;
			state.current_t = temp;
			state.chain = spawn<Chain>( engine, 
					heated_loglikelihood( temp, loglikelihood ), 
					pars_v[i%pars_v.size()], joint_prior, 1 );

			send( state.chain, atom("logger"), trace_actors[temp] );
		  chains[i] = state;
		}
	}

	void FPChainController::run() {
		size_t generation = 0;
		// Run some steps normally
		bool log_on = false;
		while (generation<warm_up+total_steps) {
			// Log to trace logger if > warm_up
			if (generation>warm_up)
				log_on = true;
			for ( auto & id_chain_state : chains ) {
				send( id_chain_state.second.chain, 
						atom("run"), no_steps_between_swaps, log_on );
			}
			// Try switching around all chains
			// Don't forget to pass correct logger around
			
			// Repeat till > total_steps + warm up
			generation += no_steps_between_swaps + 1;
		}
		
		// Close all chains
		for ( auto & id_chain_state : chains ) 
			send( id_chain_state.second.chain, atom("close") );

		size_t i = 0;
		receive_for( i, chains.size() ) (
				on( atom("closed") ) >> []() {}
		);

		// Close all loggers 
		for ( auto & temp_tr_actor : trace_actors ) 
			send( temp_tr_actor.second, atom("close") );

		i = 0;
		receive_for( i, trace_actors.size() ) (
				on( atom("closed") ) >> []() {}
		);


		// Copy traces[1] into out
		for ( auto & sample : traces[1] )
		 	out << sample << std::endl;

		// At the end gather all traces
		// Calculate means
		// Calculate marginal log likelihood
	}
};
