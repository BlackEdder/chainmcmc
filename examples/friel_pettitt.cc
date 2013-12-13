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

#include<chainmcmc/chainmcmc.hh>
#include<algorithm>
#include<cmath>

using namespace chainmcmc;

double dnorm( double x, double mu, double sigma, bool log_v = false ) {
	double pi = 3.14159265359;
	double lp = -1.0/2.0*(log(2*pi)+log(pow(sigma,2))+1.0/pow(sigma,2)*pow(x-mu,2));
	if (log_v)
		return lp;
	else
		return exp(lp);
}

double accept_prob( const std::vector<parameter_t> &new_pars,
		const std::vector<parameter_t> &old_pars, const likelihood_t &ll, 
		const std::vector<prior_t> &priors ) {
	double alpha = exp( step::mh_log_weight( ll, new_pars, priors, 1 ) -
			step::mh_log_weight( ll, old_pars, priors, 1 ) );
	//std::cout << alpha << std::endl;
	alpha = std::min( 1.0, alpha );
	//std::cout << "BLA: " << alpha << std::endl;
	return alpha;
}

double log_marginal_likelihood( std::mt19937 &eng, const likelihood_t &ll, 
		const std::vector<parameter_t> &init_pars, const 
		std::vector<prior_t> &priors ) {
	double c = 5; double n = 20;
	std::vector<double> ts;
	for (size_t i = 0; i <= n; ++i)
		ts.push_back( pow(((double) i)/n, c ) );

	double exp_ll = 0;
	double sum_ll = 0;

	double old_ti = 0;

	for ( auto & ti : ts ) {

		likelihood_t hot = [&ll, &ti]( const std::vector<parameter_t> &pars ) {
			if (ti == 0)
				return 0.0;
			else {
				return ti*ll(pars);
			}
		};
		auto chain = spawn<Chain>( eng, hot, init_pars, priors );

		std::stringstream out;
		auto logger = spawn<Logger>( out );
		send( chain, atom("logger"), logger );

		send( chain, atom("run"), 1000000, false );
		send( chain, atom("no_adapt") );
		send( chain, atom("run"), 1000000, true );
		send( chain, atom("log_weight") );
		receive( 
				on( arg_match ) >> [&out, &ti]( const double &lw ) {
				//std::cout << out.str() << std::endl;
				std::cout << ti << std::endl;
				std::cout << "Weight: " << exp( lw ) << std::endl; } );

		auto tr = trace::read_trace_per_sample( out );
		std::cout << "Prior init: " << joint_prior_t( priors )( init_pars ) << " "
			<< priors[0]( init_pars[0] ) << " " 
			<< priors[1]( init_pars[1] ) << " " 
			<< priors[2]( init_pars[2] ) << " " 
			<<std::endl;

		std::cout << tr.back() << std::endl;
		std::cout << "Prior: " << joint_prior_t( priors )( tr.back() ) << " "
			<< priors[0]( tr.back()[0] ) << " " 
			<< priors[1]( tr.back()[1] ) << " " 
			<< priors[2]( tr.back()[2] ) << " " 
			<<std::endl;

		double mean_ll = 0;
		for ( auto & sample : tr )
			mean_ll += ll( sample );
		mean_ll /= tr.size();
		std::cout << "Mean ll " << mean_ll << std::endl;

		if ( ti != 0 ) {
			sum_ll += (ti-old_ti)*(mean_ll+exp_ll)/2.0;
		}
		std::cout << "Sum ll " << sum_ll << std::endl;
		exp_ll = mean_ll;
		old_ti = ti;
	}

	return sum_ll;
}
/**
 * Goal is to implement 
 *
	Friel, N., and A. N. Pettitt. 2008. “Marginal Likelihood Estimation via Power Posteriors.” Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70 (3): 589–607. doi:10.1111/j.1467-9868.2007.00650.x.
 *
 */
int main() {
	std::mt19937 eng;
	std::vector<double> data_y = { 3040, 2470, 3610, 3480, 3810, 2330, 1800, 3110, 3160, 2310, 4360, 1880, 3670, 1740, 2250, 2650, 4970, 2620, 2900, 1670, 2540, 3840, 3800, 4600, 1900, 2530, 2920, 4990, 1670, 3310, 3450, 3600, 2850, 1590, 3770, 3850, 2480, 3570, 2620, 1890, 3030, 3030 };
	std::vector<double> data_x = { 29.2, 24.7, 32.3, 31.3, 31.5, 24.5, 19.9, 27.3, 27.1, 24.0, 33.8, 21.5, 32.2, 22.5, 27.5, 25.6, 34.5, 26.2, 26.7, 21.1, 24.1, 30.7, 32.7, 32.6, 22.1, 25.3, 30.8, 38.9, 22.1, 29.2, 30.1, 31.4, 26.7, 22.1, 30.3, 32.0, 23.2, 30.3, 29.9, 20.8, 33.2, 28.2 };
	std::vector<double> data_z = { 25.4, 22.2, 32.2, 31.0, 30.9, 23.9, 19.2, 27.2, 26.3, 23.9, 33.2, 21.0, 29.0, 22.0, 23.8, 25.3, 34.2, 25.7, 26.4, 20.0, 23.9, 30.7, 32.6, 32.5, 20.8, 23.1, 29.8, 38.1, 21.3, 28.5, 29.2, 31.4, 25.9, 21.4, 29.8, 30.6, 22.6, 30.3, 23.8, 18.4, 29.4, 28.2 };

	double mean_x = trace::details::mean_v( data_x );
	likelihood_t ll1 = [&mean_x, &data_x, &data_y]( const std::vector<parameter_t>
			&pars ) {
		double ll = 0;
		for (size_t i = 0; i < data_x.size(); ++i) {
			double mu = pars[0]+pars[1]*(data_x[i]-mean_x);
			ll += dnorm( data_y[i], mu, sqrt(pars[2]), true );
		}
		return ll;
	};

	// Second model:
	double mean_z = trace::details::mean_v( data_z );
	likelihood_t ll2 = [&mean_z, &data_z, &data_y]( const std::vector<parameter_t>
			&pars ) {
		double ll = 0;
		for (size_t i = 0; i < data_y.size(); ++i) {
			double mu = pars[0]+pars[1]*(data_z[i]-mean_z);
			ll += dnorm( data_y[i], mu, sqrt(pars[2]), true );
		}
		return ll;
	};


	std::vector<parameter_t> init_pars1 = { 3000, 185, pow(300,2) };
	std::vector<prior_t> priors1 = { prior::normal( 3000, 1e6 ),
		prior::normal( 185, 1e4 ), prior::inverse_gamma( 3.0, (2*pow(300,2)) ) };
	std::cout << "Prior init: " << joint_prior_t( priors1 )( init_pars1 ) << std::endl;

	double lml1 = log_marginal_likelihood( eng, ll1, init_pars1, priors1 );
	double lml2 = log_marginal_likelihood( eng, ll2, init_pars1, priors1 );

	std::cout << lml1 << std::endl;
	std::cout << lml2 << std::endl;
	std::cout << lml2-lml1 << " " << lml1-lml2 << std::endl;
	std::cout << exp(lml2-lml1) << " " << exp(lml1-lml2) << std::endl;

	return 0;
}
