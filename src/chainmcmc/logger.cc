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

#include "chainmcmc/logger.hh"

namespace chainmcmc {
Logger::Logger( std::ostream & out ) : out( out ) {};

void Logger::init() {
	become(
		on(atom("append"), arg_match ) >> [this]( const 
			std::string &str ) {
			out << str;
			out << std::endl;
		},
		on(atom("close")) >> []() {
			self->quit();
			return atom("closed");
		}
	);
};

TraceLogger::TraceLogger( std::vector<trace::sample_t> & tr,
		std::ostream & out ) : 
	the_trace( tr ), out( out ), log_to_stream( false ) {};

TraceLogger::TraceLogger( std::vector<trace::sample_t> & tr, bool log_to_stream,
		std::ostream &out ) :
	the_trace( tr ), out( out ), log_to_stream( log_to_stream ) {};


void TraceLogger::init() {
	become(
		on(atom("append"), arg_match ) >> [this]( const 
			std::string &str ) {
			if (log_to_stream) {
				out << str;
				out << std::endl;
			}
			trace::sample_t sample = trace::sample_from_string( str );
			the_trace.push_back( sample );
		},
		on(atom("ll"), arg_match ) >> [this]( const 
			double &ll ) {
			sum_ll += ll;
			++count_ll;
		},
		on(atom("mean_ll")) >> [this]() {
			return sum_ll/count_ll;
		},
		on(atom("close")) >> []() {
			self->quit();
			return atom("closed");
		}
	);
};
};
