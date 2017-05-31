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

#include "io.h"

#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace growing_balls {

// BEGIN class PointOfInterest
IO::PointOfInterest::PointOfInterest(std::string input_str)
    : m_elim_t(0), m_elim_partner() {
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

std::string IO::PointOfInterest::print() const {
  std::stringstream ss_result;

  ss_result.precision(22);
  ss_result << m_lat << " " << m_lon << " " << m_priority << " " << m_osmid
            << " '" << m_label << "' " << m_font_fac;

  return ss_result.str();
}

void IO::PointOfInterest::set_elimination(ElimTime elim_t, OsmId elim_p) {
  m_elim_t = elim_t;
  m_elim_partner = elim_p;
}
// END class PointOfInterest

// BEGIN class IO

bool IO::export_eliminationorder(std::string &export_file,
                                 std::vector<IO::PointOfInterest> &pois) {
  std::ofstream outFile(export_file);
  if (!outFile) {
    std::cerr << "File " << export_file << " could not be opened!" << std::endl;
    std::cerr << "CAUTION! No output file has been generated!" << std::endl;
    return false;
  }

  outFile << pois.size() << std::endl;

  outFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
  for (auto &poi : pois) {
    outFile << poi.get_lat() << ' ' << poi.get_lon() << ' ' << poi.get_osm_id()
            << ' ' << poi.get_priority() << ' ' << poi.get_elim_time() << ' '
            << poi.get_radius() << ' ' << poi.get_font_factor() << ' '
            << poi.get_label() << std::endl;
  }
  return true;
}

std::vector<IO::PointOfInterest> IO::import_label(std::string input_file) {
  std::vector<PointOfInterest> result;
  std::ifstream inFile(input_file);
  if (!inFile) {
    // if the file could not be opened: return an empty result
    std::cerr << "File " << input_file
              << " could not be opened to import label data" << std::endl;
    return result;
  }

  std::size_t total;
  inFile >> total;
  for (std::string line; std::getline(inFile, line);) {
    if (line.size() > 0 && line.at(0) != '#') {
      result.emplace_back(line);
    }
  }

  assert(result.size() == total);

  return result;
}
}

// END class IO
