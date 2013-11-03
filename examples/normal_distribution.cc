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

/**
 * Simple example where we fit a normal distribution to generated data.
 * The normal distribution has two parameters (mean and sd) that need to be fitted
 */

#include "chainmcmc/chainmcmc.hh"

using namespace chainmcmc;

std::mt19937 eng;

std::vector<double> generate_fake_data( size_t n, double mean, double sd ) {
	std::vector<double> the_data;
	std::normal_distribution<double> rnorm( mean, sd );
	for ( size_t i = 0; i<n; ++i ) {
		the_data.push_back( rnorm( eng ) );
	}
	return the_data;
}
int main() {
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

	auto contr = ChainController( loglikelihood, pars,
		priors, 10000, 100000, 4 );

	return 0;
}
