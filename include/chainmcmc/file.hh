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

#ifndef FILE_HELPERS_HH
#define FILE_HELPERS_HH

#include<vector>
#include<string>

/**
 * \brief Various functions for reading files
 */
namespace file_helpers {
	/**
	 * \brief Move to the beginning of the last n lines
	 *
	 * Implementation detail: cannot return an ifstream, not sure if move is still not
	 * implemented in c++11 or if I am missing how exactly to do it.
	 */
	void move_last_lines( std::ifstream& infile, const size_t n );

	/**
	 * \brief Get lines starting from current position and move position to the end
	 */
	std::vector<std::string> get_lines_and_move( std::ifstream& infile );

	/**
	 * \brief Read last lines from the given file
	 */
	std::vector<std::string> tail( const std::string & fname, 
			const size_t no );
};
#endif


