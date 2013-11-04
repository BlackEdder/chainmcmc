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
namespace step {
	double mh_log_weight( const likelihood_t &ll,
			const std::vector<parameter_t> &pars,
			const std::vector<prior_t> &priors, const double &temperature ) {
		double sum = 0;
		for ( size_t i = 0; i < pars.size(); ++i ) {
			double pr = priors[i]( pars[i] );
			if ( pr == 0 )
				return log(0);
			sum += log( pr );
		}
		sum += ll( pars );
		return sum/temperature;
	}

	double mh_log_weight( const double &log_likelihood,
			const std::vector<parameter_t> &pars,
			const std::vector<prior_t> &priors, 
			const double &temperature ) {
		double sum = log_likelihood;
		if (!std::isfinite( sum ))
			return sum;
		for ( size_t i = 0; i < pars.size(); ++i ) {
			double pr = priors[i]( pars[i] );
			if ( pr == 0 )
				return log(0);
			sum += log( pr );
		}
		return sum/temperature;
	}

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
			const std::vector<prior_t> &priors, double temperature ) {

		double old_lweight = mh_log_weight( old_ll, old_pars, priors, temperature );
		if (std::isnan( old_lweight ))
			return true;
		return accept( eng, old_lweight,
				mh_log_weight( ll, new_pars, priors, temperature ), temperature );
	}

	bool accept( std::mt19937 &eng, const likelihood_t &ll, 
			const std::vector<parameter_t> &old_pars,
			const std::vector<parameter_t> &new_pars,
			const std::vector<prior_t> &priors, double temperature ) {

		// First check if none of the priors are zero before we call any likelyhood
		// function
		for ( size_t i = 0; i < new_pars.size(); ++i ) {
			double pr = priors[i]( new_pars[i] );
			if ( pr == 0 )
				return false;
		}
		return accept( eng, ll, ll( old_pars ), old_pars, new_pars, 
				priors, temperature );	
	}

	parameter_t rparameter( std::mt19937 &eng, const parameter_t &par, 
			const double &sd ) {
		std::normal_distribution<double> rnorm( 0, sd );
		return par + rnorm( eng );
	}
			size_t no_accepts = 0;
			size_t no_tries = 0;
			double sd = 0.001;

	ParameterState adapt_parameter_sd( ParameterState && ps ) {
		if (ps.no_tries>100) {
			double alpha = ((double) ps.no_accepts)/ps.no_tries;
			if ( alpha > 0.25 )
				ps.sd *= 1.05;
			else if ( alpha < 0.20 )
				ps.sd /= 1.05;
			if (sd <= 0)
				ps.sd = std::numeric_limits<double>::min()*1000;
		}
		return ps;
	}

	State step( std::mt19937 &eng, State && state, const likelihood_t &ll, 
			const std::vector<prior_t> &priors, 
			double temperature ) {
		auto proposed_parameters = state.parameters;
		proposed_parameters[state.current_parameter] = 
			rparameter( eng, proposed_parameters[state.current_parameter], 
			state.pss[state.current_parameter].sd	);
		++state.pss[state.current_parameter].no_tries;

		bool accepted = false;
		double old_lweight;
		double proposed_ll;

		if (priors[state.current_parameter](
					proposed_parameters[state.current_parameter]) == 0)
			accepted = false;
		else if (std::isnan( state.loglikelihood )) {
			accepted = true;
			proposed_ll = ll( proposed_parameters );
		} else {
			proposed_ll = ll( proposed_parameters );
			old_lweight = mh_log_weight( state.loglikelihood, state.parameters,
					priors, temperature );
			if (std::isnan( old_lweight)) {
				accepted = true;
			}
			else {
				double proposed_lweight = mh_log_weight( proposed_ll, proposed_parameters,
						priors, temperature );
				accepted = accept( eng, old_lweight, proposed_lweight, temperature );
			}
		}

		if (accepted) {
			++state.pss[state.current_parameter].no_accepts;
			state.loglikelihood = proposed_ll; 
			state.parameters = proposed_parameters;
		}

		state.pss[state.current_parameter] =
			adapt_parameter_sd( std::move( state.pss[state.current_parameter] ) );
		++state.generation;
		++state.current_parameter;
		state.current_parameter = state.current_parameter%state.parameters.size();
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
	priors( priors )
{
	state.parameters = parameters;
	for ( size_t i = 0; i < parameters.size(); ++i )
		state.pss.push_back( step::ParameterState() );
}

void Chain::init()  {
	become(
		on( atom( "run" ), arg_match ) >> [this]( const int no ) {
			for (size_t i = 0; i < (size_t)no; ++i) {
				step();
			}
		},
		on( atom( "run" ), arg_match ) >> [this]( const int no, const bool start_logging ) {
			log_on = start_logging;
			for (size_t i = 0; i < (size_t)no; ++i) {
				step();
			}
		},
		on( atom( "step" ) ) >> [this]() {
			step();
		},
		on( atom("log_weight") ) >> [this]() {
			// Returns the cold weight
			return step::mh_log_weight( state.loglikelihood, state.parameters,
				priors, 1 );
		},
		on( atom("temp" ), arg_match ) >> [this]( const double &new_temp ) {
			temperature = new_temp;
		},
		on( atom("temp" ) ) >> [this]() {
			return temperature;
		},
		on( atom("close" ) ) >> []() {
			self->quit();
			return atom("closed");
		}
	);
}

void Chain::step()  {
	state = step::step( rnd_engine, std::move( state ), loglikelihood, priors,
			temperature );
	if (temperature == 1 && state.generation%50 == 0 && log_on) {
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
		s << std::endl;
		std::cout << s.str() << std::flush; 
	}
}

ChainController::ChainController( const likelihood_t &loglikelihood, 
		const std::vector<parameter_t> &parameters,
		const std::vector<prior_t> &priors, size_t warm_up, size_t total_steps,
			size_t no_chains ) : no_chains( no_chains ), warm_up( warm_up ) {

	setup( loglikelihood, {parameters}, priors, no_chains );

	run( total_steps );
}


ChainController::ChainController( const likelihood_t &loglikelihood, 
		const std::vector<std::vector<parameter_t> > &pars_v,
		const std::vector<prior_t> &priors, size_t warm_up, size_t total_steps,
		size_t no_chains ) : no_chains( no_chains ), warm_up( warm_up ) {
	setup( loglikelihood, pars_v, priors, no_chains );

	run( total_steps );
}

void ChainController::step() {
	warm_up -= 100;
	bool log_on = false;
	if (warm_up<0)
		log_on = true;
	if (chains.size()>1) {
		fisherYatesKSubsets( ids, 2, engine );
		for ( size_t i = 2; i < ids.size(); ++i ) {
			send( chains[ids[i]], atom("run"), 100, log_on );
		}

		double weight1, weight2;
		double temp1, temp2;

		size_t done = 0;
		send( chains[ids[0]], atom("log_weight") );
		receive( 
			on(arg_match) >> [&weight1,&done]( const double &ans ) {
				weight1 = ans;
				++done;
			}
		);
		send( chains[ids[1]], atom("log_weight") );
		receive( 
			on(arg_match) >> [&weight2, &done]( const double &ans ) {
				weight2 = ans;
				++done;
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
		double log_ratio = weight2/temp1+weight1/temp2 -
				weight1/temp1 - weight2/temp2;
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
			if (no_tries>10) {
				double alpha = ((double) no_accepts)/no_tries;
				if (alpha < 0.2)
					dt *= 1.05;
				if (alpha > 0.6)
					dt /= 1.05;
				if (dt == 0)
					dt = 0.01;
			}
			send( chains[ids[0]], atom("temp"), (1+dt*ids[1]) );
			send( chains[ids[1]], atom("temp"), (1+dt*ids[0]) );
		}
		send( chains[ids[1]], atom("run"), 100, log_on );
	}
	send( chains[ids[0]], atom("run"), 100, log_on );
}

	void ChainController::setup( const likelihood_t &loglikelihood, 
			const std::vector<std::vector<parameter_t> > &pars_v,
			const std::vector<prior_t> &priors,
			size_t no_chains ) {
		for ( size_t i = 0; i < no_chains; ++i ) {
			std::mt19937 eng;
			eng.seed( engine() );
			ids.push_back( i );
			double temp = (1+dt*i);
			auto chain = spawn<Chain>( eng, loglikelihood, pars_v[i%pars_v.size()],
					priors, temp );
			send( chain, atom("run"), 100 );
			chains[i] = chain;
		} 
	}
	
void ChainController::run( const size_t total_steps ) {
	for ( size_t i = 0; i < (warm_up+total_steps)/100; ++i )
		step();
	for ( auto & id_chain : chains ) 
		send( id_chain.second, atom("close") );
	size_t i = 0; 
	receive_for( i, chains.size() ) (
		on( atom("closed") ) >> []() {}
	);
}
};
