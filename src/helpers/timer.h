/*
 * Copyright 2015 Filip Krumpe <filip.krumpe@fmi.uni-stuttgart.de>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <ctime>
#include <vector>

namespace debug_timer {

class Timer
{
private:
  // ---- private class data ----
  std::chrono::time_point<std::chrono::system_clock> mStart;
  std::chrono::time_point<std::chrono::system_clock> mEnd;

  std::vector<std::chrono::time_point<std::chrono::system_clock>> mTimepoints;

public:
  // ---- constructors and destructors ----
  Timer();

  // ---- public functions ----
  void createTimepoint();

  std::vector<double> getTimes();

  double getTotal();

  void start();

  void stop();

private:
  // ---- private functions ----
  double computeDifference(
    std::chrono::time_point<std::chrono::system_clock>& t1,
    std::chrono::time_point<std::chrono::system_clock>& t2);
};
}

#endif // TIMER_H
