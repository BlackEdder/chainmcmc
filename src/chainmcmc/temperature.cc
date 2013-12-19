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

		void swap( ChainState &chainState1, ChainState &chainState2 ) {
			std::swap( chainState1.chain, chainState2.chain );
			// Send new temp
			send( chainState1.chain, atom("temp"), chainState1.current_t );
			// Send new logger
			send( chainState1.chain, atom("logger"), chainState1.logger );
			send( chainState2.chain, atom("temp"), chainState2.current_t );
			// Send new logger
			send( chainState2.chain, atom("logger"), chainState2.logger );
		}
	};
};
