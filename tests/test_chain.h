/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of chainmcmc.
  
  chainmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  chainmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with chainmcmc. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#include <cxxtest/TestSuite.h>
#include<chainmcmc/chain.hh>

using namespace chainmcmc;
using namespace chainmcmc::step;


class TestChain : public CxxTest::TestSuite 
{
	public:
		void testTemperature() {
			actor_ptr actor = spawn<Chain>();
			send( actor, atom("temp" ) );
			receive(
				on(arg_match) >> []( const double &temp ) {
					TS_ASSERT_EQUALS( temp, 1 );
				}
			);

			send( actor, atom("temp" ), 1.1 );
			send( actor, atom("temp" ) );
			receive(
				on(arg_match) >> []( const double &temp ) {
					TS_ASSERT_EQUALS( temp, 1.1 );
				}
			);
		}

		void testRun() {
			size_t count = 0;
			std::vector<parameter_t> pars = {1.0};
			std::vector<prior_t> priors = {[](const parameter_t &par){return 1;}};
			std::mt19937 rnd_eng;
			actor_ptr actor = spawn<Chain>( rnd_eng, 
					[&count]( const std::vector<parameter_t> &params ) { ++count; return 1.0; }, pars, priors );
			send( actor, atom("run"), 100 );
			send( actor, atom("temp" ) );
			receive( // Just a trick to make sure we are done with the run
				on(arg_match) >> []( const double &temp ) {
					TS_ASSERT_EQUALS( temp, 1 );
				}
			);
			TS_ASSERT_EQUALS( count, 100 );
		}

		void testAcceptance() {
			std::mt19937 eng1, eng2;
			size_t count = 0;
			likelihood_t ll = [&count]( const std::vector<parameter_t> &pars ) {
				++count;
				return pars[0];
			};
			std::vector<parameter_t> old_pars = { 1.0 };
			std::vector<parameter_t> new_pars = { 1.1 };
			prior_t prior = []( const parameter_t &par ) {
				return 0;
			};
			std::cout << "Here" << std::endl;
			TS_ASSERT( 
					!accept( eng2, ll, old_pars, new_pars, { []( const parameter_t &par ) {
				return 0;
			} } ) ); 
			TS_ASSERT_EQUALS( count, 0 ); // Make sure likelihood is not called when prior is zero
			std::cout << "Here 2" << std::endl;
			prior = []( const parameter_t &par ) {
				return 1;
			};
			TS_ASSERT( accept( eng2, ll, old_pars, new_pars,
						{ prior } ) );
			count = 0;

			TS_ASSERT( !accept( eng2, ll, old_pars, new_pars, 
			{ []( const parameter_t &par ) {
				if (par == 1.1) 
					return 0;
				else return 1;
			} } ) );
			TS_ASSERT_EQUALS( count, 0 );
			
			eng1.seed( 2 ); eng2.seed( 2 );
			size_t true_count = 0;
			for ( size_t i = 0; i<10; ++i ) {
				double ratio = 1.0/1.1;
				bool result = accept( eng2, ll, new_pars, old_pars, 
						{ prior, prior } );
				if (result)
					++true_count;
				std::uniform_real_distribution<double> unif(0,1);
				TS_ASSERT_EQUALS( result, unif(eng1) < ratio );
			}
			TS_ASSERT_EQUALS( true_count, 8 );

			// Check temperature taken into account
			eng1.seed( 2 ); eng2.seed( 2 );
			true_count = 0;
			for ( size_t i = 0; i<10; ++i ) {
				double ratio = pow(1.0/1.1, 1/1.5 );
				bool result = accept( eng2, ll, new_pars, old_pars, 
						{ prior }, 1.5 );
				if (result)
					++true_count;
				std::uniform_real_distribution<double> unif(0,1);
				TS_ASSERT_EQUALS( result, unif(eng1) < ratio );
			}
			TS_ASSERT_EQUALS( true_count, 9 );

			// Check that nan, inf and -inf are handled ok.
			TS_ASSERT( std::isnan( NAN ) );
			TS_ASSERT( 
				accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return NAN;
					else
						return 1.0;
				}, old_pars, new_pars, { prior } ) );
			// Never accept if the second solution is nan
			TS_ASSERT( 
				!accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return NAN;
					else
						return log(0);
				}, new_pars, old_pars, { prior } ) );

			TS_ASSERT( 
				!accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return NAN;
					else
						return -log(0);
				}, new_pars, old_pars, { prior } ) );

			TS_ASSERT( 
				accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return log(0);
					else
						return 1.0;
				}, old_pars, new_pars, { prior } ) );
			TS_ASSERT( 
				!accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return -log(0);
					else
						return 1.0;
				}, old_pars, new_pars, { prior } ) );
			TS_ASSERT( 
				!accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return 1.0;
					else
						return log(0);
				}, old_pars, new_pars, { prior } ) );
			TS_ASSERT( 
				accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return 1.0;
					else
						return -log(0);
				}, old_pars, new_pars, { prior } ) );
			TS_ASSERT( 
				!accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return -log(0);
					else
						return log(0);
				}, old_pars, new_pars, { prior } ) );
			TS_ASSERT( 
				accept( eng2, []( const std::vector<parameter_t> & pars ) -> double {
					if (pars[0] == 1.0) 
						return log(0);
					else
						return -log(0);
				}, old_pars, new_pars, { prior }, 1.2 ) );
		}

		void testRandomParameter() {
			std::mt19937 eng1, eng2;
			eng1.seed( 2 ); eng2.seed( 2 );
			std::normal_distribution<double> rnorm( 0, 0.001 );
			for (size_t i = 0; i < 5; ++i) {
				// This should be exactly the same, because passing engine, but is not, so using delta
				TS_ASSERT_DELTA( rparameter( eng1, 1.2, 0.001 ),
					1.2 + rnorm( eng2 ), 0.005 );
			}

			// Special rules for infinite and -infinite
			TS_ASSERT( std::isfinite( rparameter( eng1, log(0), 0.001 ) ) );
			TS_ASSERT( rparameter( eng1, log(0), 0.001 )<0 );
			TS_ASSERT( std::isfinite( rparameter( eng1, -log(0), 0.001 ) ) );
			TS_ASSERT( rparameter( eng1, -log(0), 0.001 )>0 );
		}

		void testStep() {
			// step function needs current state and returns new state.
			// Current state includes current ll, current generation, 
			// parameter index to try and current pars
			State state;
			state.parameters = { 0.5, 0.1 };
			state.pss = { ParameterState(), ParameterState() };
			// Besides that the step function needs loglikelihood function,
			// priors, temperature, rnd engine
			std::mt19937 eng1, eng2;
			eng1.seed( 2 ); eng2.seed( 2 );
			size_t count = 0;
			likelihood_t ll = [&count]( const std::vector<parameter_t> &pars ) {
				++count;
				return pars[0];
			};
			prior_t	prior = []( const parameter_t &par ) {
				return 1;
			};
			state = chainmcmc::step::step( eng2, std::move( state ), ll, 
					{ prior, prior }, true, 1.0 );
			// Test for:
			// 	generation increased, current_parameter increased
			TS_ASSERT_EQUALS( state.generation, 1 );
			TS_ASSERT_EQUALS( state.current_parameter, 1 );
			TS_ASSERT( !std::isnan(state.loglikelihood) );

			// Test parameters loop back
			state = chainmcmc::step::step( eng2, std::move( state ), ll, 
					{ prior, prior }, true, 1.0 );
			TS_ASSERT_EQUALS( state.current_parameter, 0 );
		}

		void testAdaptStepSize() {
			/*double adapt_step_size( const double current_step_size,
			const size_t no_tries, const size_t no_accepts,
			const size_t minimum_tries, const double min_level, const double max_level )*/
			TS_ASSERT_EQUALS( adapt_step_size( 1.0, 10, 5, 9, 0.3, 0.4 ),
					1.0*1.05 );
			TS_ASSERT_EQUALS( adapt_step_size( 1.0, 10, 5, 11, 0.3, 0.4 ),
					1.0 );

			TS_ASSERT_EQUALS( adapt_step_size( 1.0, 10, 2, 9, 0.3, 0.4 ),
					1.0/1.05 );

			TS_ASSERT( adapt_step_size( 0, 10, 2, 9, 0.3, 0.4 ) >	0 );
		}
};
	
