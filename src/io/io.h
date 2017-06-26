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

#ifndef IO_H
#define IO_H

#include <stdint.h>
#include <string>
#include <vector>

#include <pointofinterest.h>

namespace growing_balls {

class IO
{

public:
  static std::vector<PointOfInterest> import_label(std::string input_file);
  static bool export_eliminationorder(std::string& export_file,
                                      std::vector<PointOfInterest>& pois);
};
}

#endif // IO
