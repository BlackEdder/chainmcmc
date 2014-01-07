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
#ifndef LOGGER_HH 
#define LOGGER_HH
#include<fstream>
#include <cppa/cppa.hpp>

#include "chainmcmc/trace.hh"

namespace chainmcmc {
	using namespace cppa;

	class Logger : public event_based_actor {
		public:
			Logger( std::ostream & out );

			void init();
		protected:
			std::ostream & out;
	};

	class TraceLogger : public event_based_actor {
		public:
			/**
			 * \brief Log to trace
			 */
			TraceLogger( std::vector<trace::sample_t> & tr,
					std::ostream &out = std::cout );

			/**
			 * \brief Log to trace and also to a stream
			 */
			TraceLogger( std::vector<trace::sample_t> & tr, bool log_to_stream,
					std::ostream &out = std::cout );

			void init();
		protected:
			std::vector<trace::sample_t> & the_trace;
			std::ostream & out;
			bool log_to_stream = false;

			double sum_ll = 0;
			size_t count_ll = 0;
	};
};

#endif
