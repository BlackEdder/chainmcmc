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
#include<vector>
#include<functional>
#include<cmath>
#include<iostream>

namespace chainmcmc {
	typedef double parameter_t;
	typedef std::function<double( const parameter_t )> prior_t;

	/**
	 * \brief Joint prior, that returns prior probability for a vector of parameters
	 *
	 * Internally a vector of priors is always converted to a joint prior
	 * taking the parameters
	 *
	 * Is defined as separate class so we can easily implicetly convert between
	 * vectors of priors and a joint_prior
	 */
	class joint_prior_t {
		public:
			joint_prior_t();
			joint_prior_t( const std::function<double( const std::vector<parameter_t> 
						&pars )> &func );
			joint_prior_t( const std::vector<prior_t> &priors );
			joint_prior_t( const std::initializer_list<prior_t> &args );
			double operator()( const std::vector<parameter_t> &pars ) const;

		private:
			std::function<double( const std::vector<parameter_t> 
					&pars )> joint_func;
	};

	//typedef std::function<double( const std::vector<parameter_t> 
	//		&pars )> joint_prior_t;

	namespace prior {
		constexpr double pi() { return std::atan2(0,-1); }
		
		/**
		 * \brief Prior for normal distribution given mean and standard deviation
		 *
		 * \param sigma This is the standard deviation, not the variance
		 */
		prior_t normal( const double &mu, const double &sigma );

		prior_t inverse_gamma( const double &alpha, const double &beta );

		prior_t uniform( const double &min, const double &max );
		/**
		 * \brief Probability density function of a dirichlet distribution
		 */
		joint_prior_t dirichlet( const std::vector<double> &alpha );

};
};
#endif
