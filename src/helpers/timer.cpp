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

#include "timer.h"

debug_timer::Timer::Timer()
{
}

void
debug_timer::Timer::createTimepoint()
{
  mTimepoints.push_back(std::chrono::system_clock::now());
}

std::vector<double>
debug_timer::Timer::getTimes()
{
  std::vector<double> result;

  if (mTimepoints.size() == 0) {
    result.push_back(computeDifference(mStart, mEnd));
  } else {
    auto tp = mTimepoints.begin(), tpLast = mTimepoints.begin(),
         end = mTimepoints.end();
    result.push_back(computeDifference(mStart, *tp));

    while (std::next<>(tp) != end) {
      tpLast = tp;
      tp = ++tp;
      result.push_back(computeDifference(*tpLast, *tp));
    }
    result.push_back(computeDifference(*tp, mEnd));
  }

  return result;
}

double
debug_timer::Timer::getTotal()
{
  return computeDifference(mStart, mEnd);
}

void
debug_timer::Timer::start()
{
  mStart = std::chrono::system_clock::now();
}

void
debug_timer::Timer::stop()
{
  mEnd = std::chrono::system_clock::now();
}

double
debug_timer::Timer::computeDifference(
  std::chrono::time_point<std::chrono::_V2::system_clock>& t1,
  std::chrono::time_point<std::chrono::_V2::system_clock>& t2)
{
  std::chrono::duration<double, std::milli> difference = t2 - t1;

  return difference.count();
}
