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

#include "chainmcmc/file.hh"
#include<fstream>
#include<iostream>

namespace file_helpers {
	void move_last_lines( std::ifstream& fin, const size_t n ) {
		fin.seekg( 0,std::ios_base::end);
		if (n!=0) { 
			size_t no_found = 0;
			while (no_found < n) {
				fin.unget(); fin.unget(); // Move backwards two places

				char ch;
				fin.get(ch); // Get current char and move forward a place
				if ( (int)fin.tellg() <=1) { // We're at beginning of file
					fin.seekg(0);
					no_found = n; // Interrupt the loop
				} else {
					if (ch == '\n') {
						++no_found;
					}	
				}
			}
		}
	}

	std::vector<std::string> get_lines_and_move( std::ifstream& infile ) {
		std::vector<std::string> lines;
		std::string line;
		while (std::getline(infile, line)) {
			lines.push_back( line );
		}
		infile.clear(); // unset eofbit (and/or failbit)
		return lines;
	}

	std::vector<std::string> tail( const std::string & fname, 
			const size_t n ) {
		std::vector<std::string> lines;
		std::ifstream fin( fname );

		// Special rule that if n==0 we return whole file
		if (n!=0) {
			move_last_lines( fin, n );
		}
		
		return get_lines_and_move( fin );
	}
};

