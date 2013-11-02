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

#include<cxxtest/TestSuite.h>
#include<chainmcmc/file.hh>

#include<iostream>
#include<fstream>

using namespace file_helpers;

class TestFile : public CxxTest::TestSuite 
{
	public:
		std::string test_fname = "/tmp/file_helpers_test.txt";
		void setUp() {
			std::ofstream outFile( test_fname );
			for (size_t i = 0; i<10; ++i)
				outFile << i << std::endl;
		}

		void testMove() {
			std::ifstream infile( test_fname );
			move_last_lines( infile, 1 );
			TS_ASSERT_EQUALS( infile.tellg(), 18 );
			move_last_lines( infile, 2 );
			TS_ASSERT_EQUALS( infile.tellg(), 16 );

			move_last_lines( infile, 100 );
			TS_ASSERT_EQUALS( infile.tellg(), 0 );
			move_last_lines( infile, 0 );
			TS_ASSERT_EQUALS( infile.tellg(), 20 );
		}

		void testGetAndMove() {
			std::ifstream infile( test_fname );
			move_last_lines( infile, 2 );
			auto lines = get_lines_and_move( infile );
			TS_ASSERT_EQUALS( lines.size(), 2 );
			TS_ASSERT_EQUALS( infile.tellg(), 20 );


			lines = get_lines_and_move( infile );
			TS_ASSERT_EQUALS( lines.size(), 0 );

			std::ofstream outFile( test_fname, std::fstream::app );
			for (size_t i = 0; i<10; ++i)
				outFile << i << std::endl;

			lines = get_lines_and_move( infile );
			TS_ASSERT_EQUALS( lines.size(), 10 );
		}

		void testTail() {
			auto lines = tail( test_fname, 1 );
			TS_ASSERT_EQUALS( lines.size(), 1 );
			TS_ASSERT_EQUALS( lines[0], "9" );
			lines = tail( test_fname, 2 );
			TS_ASSERT_EQUALS( lines.size(), 2 );
			TS_ASSERT_EQUALS( lines[0], "8" );
			TS_ASSERT_EQUALS( lines[1], "9" );

			lines = tail( test_fname, 0 );
			TS_ASSERT_EQUALS( lines.size(), 10 );

			lines = tail( test_fname, 100 );
			TS_ASSERT_EQUALS( lines.size(), 10 );
		}
};
