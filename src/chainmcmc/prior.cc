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

#include<cassert>

#include "chainmcmc/prior.hh"

#include <boost/math/distributions/inverse_gamma.hpp>

namespace chainmcmc {

	/*
	 * joint prior class
	 */
	joint_prior_t::joint_prior_t() {}
	joint_prior_t::joint_prior_t( 
			const std::function<double( const std::vector<parameter_t> 	&pars )> &func ) :
		joint_func( func ) {}

	joint_prior_t::joint_prior_t( const std::vector<prior_t> &priors )
	{
		joint_func =  [&priors]( const std::vector<parameter_t> &pars ) {
			assert( pars.size() <= priors.size() );
			double prob = 1;
			for ( size_t i = 0; i < pars.size(); ++i ) {
				double pr = priors[i](pars[i]);
				if (pr == 0)
					return 0.0;
				else
					prob *= pr;
			}
			return prob;
		};
	}

	joint_prior_t::joint_prior_t( const std::initializer_list<prior_t> &args ) {
		joint_func =  [&args]( const std::vector<parameter_t> &pars ) {
			assert( pars.size() <= args.size() );

			double prob = 1;
			auto it = args.begin();
			for ( size_t i = 0; i < pars.size(); ++i ) {
				double pr = (*it)(pars[i]);
				if (pr == 0)
					return 0.0;
				else
					prob *= pr;
				++it;
			}
			return prob;
		};
	}

	double joint_prior_t::operator()( const std::vector<parameter_t> &pars ) const {
		return joint_func( pars );
	}

	namespace prior {

		prior_t normal( const double &mu, const double &sigma ) {
			return [mu, sigma]( const parameter_t &x ) {
				return 1.0/(sigma*sqrt(2.0*pi()))*exp(-pow(x-mu,2)/(2*pow(sigma,2)));
			};
		}

		prior_t inverse_gamma( const double &alpha, const double &beta ) {
			return [alpha, beta]( const parameter_t &x ) {
				if (x>=0)
					return boost::math::pdf( 
						boost::math::inverse_gamma_distribution<double>( alpha,
							beta ), x );
				else
					return 0.0;
			};
		}

		prior_t uniform( const double &min, const double &max ) {
			double prob = 1.0/(max-min); // No need to calculate everytime
			return [min, max, prob]( const double &x ) {
				if (x>=min && x<=max)
					return prob;
				else
					return 0.0;
			};
		}

		joint_prior_t dirichlet( const std::vector<double> &alphas ) {
			// Calculate weight once.
			double weight = 1;
			double sum_alphas = 0;
			for ( auto & alpha : alphas ) {
				weight /= boost::math::tgamma( alpha ); 
				sum_alphas += alpha;
			}
			weight *= boost::math::tgamma( sum_alphas );


			return joint_prior_t( [alphas, weight]( 
						const std::vector<parameter_t> &pars ) {
				if (alphas.size()-1 != pars.size())
					return 0.0;
				double sum = 0;
				for ( auto &par : pars ) {
					if (par<0 || par>1)
						return 0.0;
					sum += par;
				}
				if (sum>1)
					return 0.0;

				double prob = 1.0;
				for ( size_t i = 0; i < alphas.size(); ++i ) {
					if ( i<pars.size() )
						prob *= pow(pars[i],alphas[i]-1);
					else
						prob *= pow( 1-sum, alphas[i]-1);
				}
				return prob*weight;
			} );
		}
	};
};
