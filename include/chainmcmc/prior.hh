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

#ifndef PRIOR_HH 
#define PRIOR_HH
#include<functional>
#include<cmath>

namespace chainmcmc {
	typedef double parameter_t;
	typedef std::function<double( const parameter_t )> prior_t;
	namespace prior {
		constexpr double pi() { return std::atan2(0,-1); }
		
		/**
		 * \brief Prior for normal distribution given mean and standard deviation
		 *
		 * \param sigma This is the standard deviation, not the variance
		 */
		prior_t normal( const double &mu, const double &sigma );
	};
};
#endif
