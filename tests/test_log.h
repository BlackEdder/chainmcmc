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

#include<sstream>
#include<cxxtest/TestSuite.h>
#include<chainmcmc/chainmcmc.hh>

using namespace chainmcmc;
using namespace cppa;

class TestLog : public CxxTest::TestSuite 
{
	public:
		void testAppend() {
			std::ostringstream s;
			auto logger = spawn<Logger>( s );

			send( logger, atom("append"), "1.2\t0.9" );
			send( logger, atom("close") );

			receive(
				on( atom("closed") ) >> []() {}
			);	
			TS_ASSERT_EQUALS( s.str(), "1.2\t0.9\n" );
		}
};
