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

#include "chainmcmc/temperature.hh"

namespace chainmcmc {
	namespace temperature {
		bool ChainState::operator<( const ChainState &cs ) {
			return current_t < cs.current_t;
		}

		bool accept( std::mt19937 &engine, const ChainState &chainState1, 
				const ChainState &chainState2 ) {
			// Make sure we are finished running
			double ll1; double ljp1;
			send( chainState1.chain, atom("llikelih") );
			receive( 
				on(arg_match) >> [&ll1]( const double &ll ) {
					ll1 = ll;
				}
			);
			send( chainState1.chain, atom("lprior") );
			receive( 
				on(arg_match) >> [&ljp1]( const double &ljp ) {
					ljp1 = ljp;
				}
			);

			double ll2; double ljp2;
			send( chainState2.chain, atom("llikelih") );
			receive( 
				on(arg_match) >> [&ll2]( const double &ll ) {
					ll2 = ll;
				}
			);
			send( chainState2.chain, atom("lprior") );
			receive( 
				on(arg_match) >> [&ljp2]( const double &ljp ) {
					ljp2 = ljp;
				}
			);


			double alpha = 
					chainState2.current_t * ll1 + ljp1 +
					chainState1.current_t * ll2 + ljp2 - (
						chainState1.current_t * ll1 + ljp1 +
						chainState2.current_t * ll2 + ljp2 
					);
			alpha = exp( alpha );
			bool accepted = false;
			if ( alpha > 1 )
				accepted = true;
			else {
				std::uniform_real_distribution<double> unif(0,1);
				if ( unif( engine ) < alpha )
					accepted = true;
			}
			return accepted;
		}

		void swap( ChainState &chainState1, ChainState &chainState2 ) {
			std::swap( chainState1.chain, chainState2.chain );
			// Send new temp
			send( chainState1.chain, atom("temp"), chainState1.current_t );
			send( chainState2.chain, atom("temp"), chainState2.current_t );
			// Send new logger
			send( chainState1.chain, atom("logger"), chainState1.logger );
			send( chainState2.chain, atom("logger"), chainState2.logger );
		}
	};
};
