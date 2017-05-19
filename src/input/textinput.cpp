/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  Filip Krumpe <filip.krumpe@fmi.uni-stuttgart.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "textinput.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>

// BEGIN class PointOfInterest
growing_balls::TextInput::PointOfInterest::PointOfInterest(std::string input_str)
{
  std::stringstream ss(input_str);
  std::string first, second, last; 
  std::getline(ss, first, '\'');
  std::getline(ss, second, '\'');
  std::getline(ss, last);
  
  m_label = second;
  std::stringstream ss_first(first);
  ss_first >> m_lat >> m_lon >> m_priority >> m_radius >> m_osmid;
  std::stringstream ss_last(last);
  ss_last >> m_font_fac;
}

std::string growing_balls::TextInput::PointOfInterest::print() const
{
  std::stringstream ss_result;
  
  ss_result.precision(22);
  ss_result << m_lat << " " << m_lon << " " << m_priority << " " << m_osmid << " '" << m_label << 
  "' " << m_font_fac;
  
  return ss_result.str();
}
// END class PointOfInterest

// BEGIN class TextInput
std::vector<growing_balls::TextInput::PointOfInterest> growing_balls::TextInput::import_label(std::string input_file)
{
  std::vector<PointOfInterest> result;
  std::ifstream inFile(input_file);
  if (!inFile) {
    // if the file could not be opened: return an empty result
    std::cerr << "File " << input_file << " could not be opened to import label data" << std::endl;
    return result;
  }
  
  std::size_t total;
  inFile >> total;
  for (std::string line; std::getline(inFile, line);) {
    if (line.size() > 0 && line.at(0) != '#') {
      result.emplace_back(line);
    }
  }
  
  assert (result.size() == total);
  
  return result;
}
// END class TextInput
