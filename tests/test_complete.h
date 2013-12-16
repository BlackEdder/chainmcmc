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
#include<chainmcmc/chainmcmc.hh>

using namespace chainmcmc;


class TestComplete : public CxxTest::TestSuite 
{
	public:
		std::mt19937 eng;


		std::vector<double> generate_fake_data( size_t n, double mean, double sd ) {
			std::vector<double> the_data;
			std::normal_distribution<double> rnorm( mean, sd );
			for ( size_t i = 0; i<n; ++i ) {
				the_data.push_back( rnorm( eng ) );
			}
			return the_data;
		}

		void testNormalDistribution() {
			// generate data
			auto the_data = generate_fake_data( 100, -1, 0.2 );

			std::vector<parameter_t> pars  = { 1, 1 };

			auto loglikelihood = [&the_data]( const std::vector<parameter_t> &params ) {
				double pi = 3.1415926535897;
				double mean = params[0];
				double sd = params[1];
				double n = the_data.size();
				double sum = 0;
				for ( auto &d : the_data )
					sum += pow(d-mean,2);
				double ll = -n/2.0*(log(2*pi)+log(pow(sd,2))) - 1.0/(2.0*pow(sd,2))*sum;
				return ll;
			};

			std::vector<prior_t> priors = { 
				[]( const parameter_t &par ) { return 1; }, //Any value is allowed for mean
				[]( const parameter_t &par ) { if (par<0) return 0; else return 1; } // sd need to be positive
			};

			// Use chain and log to ostringstream
			std::stringstream output;
			auto contr = ChainController( loglikelihood, pars,
					priors, 10000, 10000, 4, output );

			// Use trace to analyse ostringstream
			auto tr = trace::read_trace_per_sample( output );
			// Check results
			
			auto means = trace::means( tr );

			TS_ASSERT_DELTA( means[0], -1, 0.01 );
			TS_ASSERT_DELTA( means[1], 0.2, 0.05 );
		}

		void testFPNormalDistribution() {
			// generate data
			auto the_data = generate_fake_data( 100, -1, 0.2 );

			std::vector<parameter_t> pars  = { 1, 1 };

			auto loglikelihood = [&the_data]( const std::vector<parameter_t> &params ) {
				double pi = 3.1415926535897;
				double mean = params[0];
				double sd = params[1];
				double n = the_data.size();
				double sum = 0;
				for ( auto &d : the_data )
					sum += pow(d-mean,2);
				double ll = -n/2.0*(log(2*pi)+log(pow(sd,2))) - 1.0/(2.0*pow(sd,2))*sum;
				return ll;
			};

			std::vector<prior_t> priors = { 
				[]( const parameter_t &par ) { return 1; }, //Any value is allowed for mean
				[]( const parameter_t &par ) { if (par<0) return 0; else return 1; } // sd need to be positive
			};

			// Use chain and log to ostringstream
			std::stringstream output;
			auto contr = FPChainController( loglikelihood, { pars },
			  priors, 10000, 50000, output );
			contr.run();

			// Use trace to analyse ostringstream
			auto tr = trace::read_trace_per_sample( output );
			// Check results
			
			auto means = trace::means( tr );

			TS_ASSERT_DELTA( means[0], -1, 0.02 );
			TS_ASSERT_DELTA( means[1], 0.2, 0.05 );
		}
};


