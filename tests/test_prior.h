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

#include <cxxtest/TestSuite.h>
#include<chainmcmc/chainmcmc.hh>

using namespace chainmcmc;
using namespace chainmcmc::prior;


class TestPrior : public CxxTest::TestSuite 
{
	public:
		void testNormal() {
			auto correct = []( const double &x, const double &mu, const double &sigma ) {
				return 1.0/(sigma*sqrt(2*pi()))*exp(-pow(x-mu,2)/(2*pow(sigma,2)));
			};

			auto pr1 = prior::normal( 1, 0.1 );
			TS_ASSERT_EQUALS( correct( 1, 1, 0.1 ), pr1( 1 ) );
			TS_ASSERT_EQUALS( correct( 1, 1, 0.1 ), pr1( 0.1 ) );
			TS_ASSERT_DIFFERS( correct( 1, 1, 0.01 ), pr1( 1 ) );
			auto pr1 = prior::normal( 2, 0.01 );
			TS_ASSERT_EQUALS( correct( 1, 2, 0.01 ), pr1( 1 ) );
			TS_ASSERT_EQUALS( correct( 0.1, 2, 0.01 ), pr1( 0.1 ) );
			TS_ASSERT_DIFFERS( correct( 1, 2, 0.1 ), pr1( 1 ) );
		}
};
