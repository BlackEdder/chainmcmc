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

#ifndef TRACE_HH
#define TRACE_HH
#include<math.h>
#include<numeric>
#include<algorithm>
#include<vector>
#include<iostream>

#include <boost/lexical_cast.hpp>
#include <csv_parser/csv_parser.hpp>

#include "chainmcmc/chain.hh"

namespace chainmcmc {
namespace trace {
	typedef std::vector<parameter_t> sample_t;


	namespace details {
		double mean_v( const std::vector<double> &v ) {
			double sum = std::accumulate( v.begin(), v.end(), 0.0, 
					[]( const double a, const double b ) { return a+b; } );
			return sum/v.size();
		}

		double var_v( const std::vector<double> &v ) {
			double mean = mean_v( v );
			double var_sum = std::accumulate( v.begin(), v.end(), 0.0, 
					[&mean]( const double a, const double b ) { 
					return a+pow(b-mean,2); } );
			return var_sum/v.size();
		}

		std::pair<double, double> confidence( const double interval,
				std::vector<double> v ) {
			std::pair<double, double> result;
			std::sort( v.begin(), v.end() );
			double alpha = (1.0-interval)/2.0;
			result.first = v[ceil(alpha*v.size())];
			result.second = v[floor((1.0-alpha)*v.size())];
			return result;
		}

		double cov_v( const std::vector<double> &v,
				const std::vector<double> &w ) {
			if ( v.size() != w.size() ) {
				std::cerr << "Vectors need to be the same size" << std::endl;
				throw;
			}
			double meanv = mean_v( v );
			double meanw = mean_v( w );
			double cov = 0;
			for ( size_t i = 0; i < v.size(); ++i ) {
				cov += (v[i] - meanv)*(w[i]-meanw);
			}
			return cov/v.size();
		}
	};

	/**
	 * \brief Reads a trace file
	 *
	 * Returns a vector with vectors of each parameter value
	 */
	std::vector<std::vector<double> > read_trace( 
			const std::string & fname );

	/**
	 * \brief Reads a trace file
	 *
	 * Returns a vector with samples
	 *
	 * \param tail Read the given number of samples from the end. Set to zero to read the whole file
	 */
	std::vector<std::vector<parameter_t> > read_trace_per_sample( 
			std::istream& infile, const size_t tail = 0 );

	/**
	 * \brief Read any new samples that have been appended to the trace file
	 */
	std::vector<std::vector<parameter_t> > follow_trace_per_sample( 
			std::istream& infile );


	/**
	 * \brief Returns randomly choosed samples from the trace file
	 *
	 * \param tail Only use the tail end of the file. If 0 then whole file is used
	 */
	std::vector<std::vector<parameter_t> > random_samples_from_file( 
			const size_t &no, const std::string & filename, const size_t tail = 0 ); 

	/**
	 * \brief Returns means of the parameters in the trace
	 *
	 * Trace is assumed to be a vector of samples
	 */
	std::vector<double> means( const std::vector<sample_t> &samples );

	std::vector<double> variances_sample( 
			const std::vector<sample_t> & samples );
};
};
#endif
