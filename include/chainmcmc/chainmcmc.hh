/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Edwin van Leeuwen
  
  This file is part of libmcmc.
  
 	libmcmc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  libmcmc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Gillespie. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef MCMC_HH
#define MCMC_HH
#include "chainmcmc/parameter.hh"
#include "chainmcmc/trace.hh"
#include "chainmcmc/temperature.hh"
#include "chainmcmc/chain.hh"
#include "chainmcmc/logger.hh"
#include "chainmcmc/prior.hh"
/**
 * \mainpage MCMC library for model fitting
 *
 * If you are looking at finding the posterior distributions for different
 * parameter values then you should basically use the chainmcmc::ChainController class.
 *
 * For model comparison use the chainmcmc::FPChainController class which is slower, but
 * can approximate the marginal likelihood of your model using the method 
 * described in:
 Friel, N., and A. N. Pettitt. 2008. “Marginal Likelihood Estimation via Power Posteriors.” Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70 (3): 589–607. doi:10.1111/j.1467-9868.2007.00650.x.
 *
 * Examples are provided in the examples directory. Another good place to look
 * is tests/test_complete.h
 */


/**
 * \brief MCMC library for model fitting
 *
 * If you are looking at finding the posterior distributions for different
 * parameter values then you should basically use the ChainController class.
 *
 * For model comparison use the FPChainController class is slower, but can 
 * approximate the marginal likelihood of your model using the method 
 * described in:
 Friel, N., and A. N. Pettitt. 2008. “Marginal Likelihood Estimation via Power Posteriors.” Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70 (3): 589–607. doi:10.1111/j.1467-9868.2007.00650.x.
 */
namespace chainmcmc {
};
#endif
