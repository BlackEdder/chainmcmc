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
	alpha = std::min( 1.0, alpha );
	return alpha;
}

double marginal_likelihood( const likelihood_t &ll, 
		const std::vector<prior_t> &priors, 
		const std::vector<std::vector<parameter_t> > &tr ) {
	std::mt19937 eng;
	auto theta_star = trace::means( tr );
	auto theta_var = trace::variances_sample( tr );
	size_t no_samples = 1000;
	double numer = 0;
	double denom = 0;
	for (size_t i = 0; i < no_samples; ++i) {
		auto theta_g = tr[ tr.size()-1-i ];
		double tmp = accept_prob( theta_g, theta_star, ll, priors );
		for ( size_t j = 0; j < theta_star.size(); ++j )
			tmp *= dnorm( theta_g[j], theta_star[j], theta_var[j], false );
		numer += tmp;

		theta_g = theta_star;
		for ( size_t j = 0; j < theta_star.size(); ++j ) {
			std::normal_distribution<double> rnorm( 0, theta_var[j] );
			theta_g[j] += rnorm( eng );
		}

		denom += accept_prob( theta_star, theta_g, ll, priors );
	}
	double lpi = log( numer/denom );

	return exp(step::mh_log_weight( ll, theta_star, priors, 1 )-lpi);
}

/**
 * Goal is to implement first example from:
 * Han, Cong, and Bradley P Carlin. “Markov Chain Monte Carlo Methods for Computing Bayes Factors: A Comparative Review.” Journal of the American Statistical Association 96, no. 455 (September 2001): 1122–1132. doi:10.1198/016214501753208780.
 *
 * And test different methods for correct implementation
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

	std::vector<parameter_t> init_pars1 = { 3000, 185, pow(300,2) };
	std::vector<prior_t> priors1 = { prior::normal( 3000, 1e6 ),
		prior::normal( 185, 1e4 ), prior::inverse_gamma( 3.0, 1.0/(2*pow(300,2)) ) };

	auto chain = spawn<Chain>( eng, ll1, init_pars1, priors1 );

	std::stringstream out;
	auto logger = spawn<Logger>( out );
	send( chain, atom("logger"), logger );

	send( chain, atom("run"), 1000000, false );
	send( chain, atom("run"), 1000000, true );
	send( chain, atom("log_weight") );
	receive( 
			on( arg_match ) >> [&out]( const double &lw ) {
			//std::cout << out.str() << std::endl;
			std::cout << "Weight: " << exp( lw ) << std::endl; } );

	auto tr1 = trace::read_trace_per_sample( out );

	// Calculate marginal likelihood
	std::cout << marginal_likelihood( ll1, priors1, tr1 ) << std::endl;

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

	chain = spawn<Chain>( eng, ll2, init_pars1, priors1 );

	out.clear();
	logger = spawn<Logger>( out );
	send( chain, atom("logger"), logger );

	send( chain, atom("run"), 1000000, false );
	send( chain, atom("run"), 1000000, true );
	send( chain, atom("log_weight") );
	receive( 
			on( arg_match ) >> [&out]( const double &lw ) {
			//std::cout << out.str() << std::endl;
			std::cout << "Weight: " << exp( lw ) << std::endl; } );

	auto tr2 = trace::read_trace_per_sample( out );

	auto means_tr1 = trace::means( tr1 );
	auto means_tr2 = trace::means( tr2 );
	for (size_t i = 0; i < means_tr1.size(); ++i)
		std::cout << means_tr1[i] << " " << means_tr2[i] << std::endl;
	std::cout << marginal_likelihood( ll2, priors1, tr2 ) << std::endl;

	return 0;
}
