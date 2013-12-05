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
			TS_ASSERT_EQUALS( correct( 0.1, 1, 0.1 ), pr1( 0.1 ) );
			TS_ASSERT_DIFFERS( correct( 1, 1, 0.01 ), pr1( 1 ) );
			pr1 = prior::normal( 2, 0.01 );
			TS_ASSERT_EQUALS( correct( 1, 2, 0.01 ), pr1( 1 ) );
			TS_ASSERT_EQUALS( correct( 0.1, 2, 0.01 ), pr1( 0.1 ) );
			TS_ASSERT_DIFFERS( correct( 1, 2, 0.1 ), pr1( 1 ) );
		}

		void testIG() {
			auto pr = prior::inverse_gamma( 3, 1 );
			// Proper values are taken from wikipedia plot (not very precise)
			TS_ASSERT_DELTA( pr( 0.5 ), 1.1, 0.05 );
			TS_ASSERT_DELTA( pr( 0.25 ), 2.3, 0.05 );
			pr = prior::inverse_gamma( 1, 1 );
			TS_ASSERT_DELTA( pr( 0.5 ), 0.5, 0.05 );
			TS_ASSERT_EQUALS( pr( -0.5 ), 0.0 );
		}

		void testUnif() {
			auto pr = prior::uniform( -1, 9 );
			TS_ASSERT_EQUALS( pr( 0 ), 0.1 );
			TS_ASSERT_EQUALS( pr( 5 ), 0.1 );
			TS_ASSERT_EQUALS( pr( -2 ), 0 );
			TS_ASSERT_EQUALS( pr( 10 ), 0 );
		}

		void testDirichlet() {
			// Proper values are taken from Mathematica
			auto jp = prior::dirichlet( {2,3,2} );
			TS_ASSERT_EQUALS( 0, jp({ 0.5, 0.6 }) );
			TS_ASSERT_EQUALS( 0, jp({ -0.5, 0.6 }) );
			TS_ASSERT_EQUALS( 0, jp({ 0.5, -0.3 }) );

			TS_ASSERT_DELTA( 3.24, jp({ 0.5, 0.3 }), 0.01 );
			jp = prior::dirichlet( {2,1,2} );
			TS_ASSERT_DELTA( 2.4, jp({ 0.5, 0.3 }), 0.01 );
			TS_ASSERT_DELTA( 1.44, jp({ 0.3, 0.5 }), 0.01 );
		}
};
