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

#include "chainmcmc/prior.hh"

#include <boost/math/distributions/inverse_gamma.hpp>

namespace chainmcmc {
	namespace prior {
		prior_t normal( const double &mu, const double &sigma ) {
			return [&mu, &sigma]( const parameter_t &x ) {
				return 1.0/(sigma*sqrt(2.0*pi()))*exp(-pow(x-mu,2)/(2*pow(sigma,2)));
			};
		}

		prior_t inverse_gamma( const double &alpha, const double &beta ) {
			return [&alpha, &beta]( const parameter_t &x ) {
				if (x>=0)
					return boost::math::pdf( 
						boost::math::inverse_gamma_distribution<double>( alpha,
							beta ), x );
				else
					return 0.0;
			};
		}

		prior_t uniform( const double &min, const double &max ) {
			static double prob = 1.0/(max-min); // No need to calculate everytime
			return [&min, &max]( const double &x ) {
				if (x>=min && x<=max)
					return prob;
				else
					return 0.0;
			};
		}
	};
};