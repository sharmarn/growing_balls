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

#ifndef ELIMINATIONORDER_H
#define ELIMINATIONORDER_H

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <queue>
#include <stdint.h>
#include <string.h>
#include <unordered_map>

#include "geofunctions.h"
#include "io.h"
#include "pointofinterest.h"
#include "spatialhelper.h"
#include "timer.h"

namespace {

using ID = growing_balls::OsmId;
using PoiMap = std::unordered_map<ID, growing_balls::PointOfInterest>;

enum class Heuristic
{
  DEFAULT,
  RADIUS,
  RANDOM,
  IN_RANGE
};

// Set a heuristic from above to compute the elimination order with.
Heuristic choose_heuristic = Heuristic::IN_RANGE;

/*
 * This function is called, when the non-heuristical approach 'DEFAULT' is set.
 * Always prefer the center with the greater OSM ID.
 */
bool
prefer_center_with_greater_osm_id(const growing_balls::PointOfInterest& p1,
                                  const growing_balls::PointOfInterest& p2)
{
  if (p1.get_osm_id() > p2.get_osm_id()) {
    return true;
  }

  else
    return false;
};

/*
 * Always prefer the center with the smaller radius.
 * If the radii of both centers are equal, then choose the center with the
 * greater OSM ID.
 */
bool
prefer_center_with_smaller_radius(const growing_balls::PointOfInterest& p1,
                                  const growing_balls::PointOfInterest& p2)
{
  if (p1.get_radius() < p2.get_radius()) {
    return true;
  }

  else if (p1.get_radius() == p2.get_radius()) {
    if (p1.get_osm_id() > p2.get_osm_id()) {
      return true;
    }

    else
      return false;
  }

  else
    return false;
};

int
flip()
{
  return rand() % 2;
};

// Always choose one of the centers randomly.
bool
prefer_center_randomly()
{
  int coin = 0;
  coin = flip();
  if (coin == 0) {
    return true;
  }

  else
    return false;
};

/*
 * Always prefer the center, that has lesser centers in range 2*|c_i nn|,
 * with nn being the nearest neighbor. If both centers have the equal amount of
 * centers in range, then choose the center with the greater osm ID.
 */
bool
prefer_center_with_lesser_centers_in_range(
  const growing_balls::PointOfInterest& p1,
  const growing_balls::PointOfInterest& p2,
  growing_balls::SpatialHelper& sh,
  const PoiMap& poi_map)
{
  auto id_nn_of_p1 = sh.get_nearest_neighbor(p1.get_osm_id());
  auto& nn_p1 = poi_map.at(id_nn_of_p1);

  auto id_nn_of_p2 = sh.get_nearest_neighbor(p2.get_osm_id());
  auto& nn_p2 = poi_map.at(id_nn_of_p2);

  auto distance_p1_nn = growing_balls::distance_in_centimeters(
    p1.get_lat(), p1.get_lon(), nn_p1.get_lat(), nn_p1.get_lon());

  auto distance_p2_nn = growing_balls::distance_in_centimeters(
    p2.get_lat(), p2.get_lon(), nn_p2.get_lat(), nn_p2.get_lon());

  int count_p1_in_range_centers = 0;
  for (auto id : sh.get_in_range(p1.get_osm_id(), 2 * distance_p1_nn)) {
    auto& p1_in_range = poi_map.at(id);
    if (p1.get_priority() >= p1_in_range.get_priority()) {
      count_p1_in_range_centers++;
    }
  }

  int count_p2_in_range_centers = 0;
  for (auto id : sh.get_in_range(p2.get_osm_id(), 2 * distance_p2_nn)) {
    auto& p2_in_range = poi_map.at(id);
    if (p2.get_priority() >= p2_in_range.get_priority()) {
      count_p2_in_range_centers++;
    }
  }

  if (count_p1_in_range_centers < count_p2_in_range_centers) {
    return true;
  }

  else if (count_p1_in_range_centers == count_p2_in_range_centers) {
    if (p1.get_osm_id() > p2.get_osm_id()) {
      return true;
    }

    else
      return false;
  }

  else
    return false;
};

bool
use_heuristic(const growing_balls::PointOfInterest& p1,
              const growing_balls::PointOfInterest& p2,
              growing_balls::SpatialHelper& sh,
              const PoiMap& poi_map)
{
  assert(p1.get_priority() == p2.get_priority());

  if (choose_heuristic == Heuristic::RADIUS) {
    return prefer_center_with_smaller_radius(p1, p2);
  }

  else if (choose_heuristic == Heuristic::IN_RANGE) {
    return prefer_center_with_lesser_centers_in_range(p1, p2, sh, poi_map);
  }

  else if (choose_heuristic == Heuristic::RANDOM) {
    return prefer_center_randomly();
  }

  else
    return prefer_center_with_greater_osm_id(p1, p2);
};

bool
prefer_p1(const growing_balls::PointOfInterest& p1,
          const growing_balls::PointOfInterest& p2,
          growing_balls::SpatialHelper& sh,
          const PoiMap& poi_map)
{
  if (p1.get_priority() > p2.get_priority()) {
    return true;
  }

  else if (p1.get_priority() == p2.get_priority()) {
    return use_heuristic(p1, p2, sh, poi_map);
  }

  else
    return false;
};
}

namespace growing_balls {
using Time = ElimTime;

class EliminationOrder
{
public:
  struct Elimination
  {
    PointOfInterest m_eliminated;
    PointOfInterest m_eliminated_by;
    Time m_elimination_time;

    Elimination(Time elim_t, PointOfInterest elim, PointOfInterest elim_by)
      : m_eliminated(std::move(elim))
      , m_eliminated_by(std::move(elim_by))
      , m_elimination_time(elim_t){};
  };

public:
  EliminationOrder() = default;

  EliminationOrder(EliminationOrder&& other) = default;
  EliminationOrder& operator=(EliminationOrder&& other) = default;

  std::vector<Elimination> compute_elimination_order(std::string file);

private:
  std::vector<double> m_debug_times;
};
}

#endif // ELIMINATIONORDER_H

// BEGIN DEFINITION

namespace growing_balls {

// BEGIN Helpers
namespace {
using ID = OsmId;
using PoiMap = std::unordered_map<ID, PointOfInterest>;

enum class EventType : int32_t
{
  UPDATE_EVENT = 1,
  COLLISION_EVENT = 2,
};

struct Event
{
  Time m_trigger_time;
  ID m_coll1;
  ID m_coll2;

  EventType m_evt_type;

  Event(Time t, ID c1)
    : m_trigger_time(t)
    , m_coll1(c1)
    , m_coll2(c1)
    , m_evt_type(EventType::UPDATE_EVENT){};

  Event(Time t, ID c1, ID c2)
    : m_trigger_time(t)
    , m_coll1(c1)
    , m_coll2(c2)
    , m_evt_type(EventType::COLLISION_EVENT){};

  // In order to use a priority_queue (a max heap) as min heap: overload
  // operator <
  bool operator<(const Event& other) const
  {
    if (m_trigger_time == other.m_trigger_time) {
      if (m_evt_type == EventType::UPDATE_EVENT) {
        // update events are prefered to collision events
        return false;
      } else if (m_evt_type == other.m_evt_type) {
        // make things deterministic: use m_coll1 ids to break ties with equal
        // event type
        return m_coll1 > other.m_coll1;
      } else {
        // other is an update event and we are not -> other has higher priority
        return true;
      }
    } else {
      return m_trigger_time > other.m_trigger_time;
    }
  }
};

Time
compute_collision_time(const PointOfInterest& p1, const PointOfInterest& p2)
{
  auto d = distance_in_centimeters(
    p1.get_lat(), p1.get_lon(), p2.get_lat(), p2.get_lon());

  return d / (p1.get_radius() + p2.get_radius());
}

Event
predict_collision(const PointOfInterest& p,
                  Time t,
                  SpatialHelper& sh,
                  const PoiMap& poi_map)
{
  auto id_nn = sh.get_nearest_neighbor(p.get_osm_id());

  if (id_nn == sh.UNDEFINED_ID) {
    return Event(
      std::numeric_limits<Time>::max(), p.get_osm_id(), p.get_osm_id());
  }

  auto& nn = poi_map.at(id_nn);

  auto distance_pnn = distance_in_centimeters(
    p.get_lat(), p.get_lon(), nn.get_lat(), nn.get_lon());
  Time t_upd = distance_pnn / (2 * p.get_radius());

  if (t < t_upd) {
    return Event(t_upd, p.get_osm_id());
  } else {
    Time min_coll_t = std::numeric_limits<Time>::max();
    OsmId min_coll_p = 0;
    for (auto id : sh.get_in_range(p.get_osm_id(), 2 * distance_pnn)) {
      auto& p2 = poi_map.at(id);

      Time coll_t = compute_collision_time(p, p2);
      if (coll_t < min_coll_t) {
        min_coll_t = coll_t;
        min_coll_p = p2.get_osm_id();
      }
    }

    return Event(min_coll_t, p.get_osm_id(), min_coll_p);
  }
}
}
// END Helpers

// BEGIN class EliminationOrder

std::vector<EliminationOrder::Elimination>
EliminationOrder::compute_elimination_order(std::string file)
{
  debug_timer::Timer timer;
  timer.start();
  auto labels = IO::import_label(file);

  timer.createTimepoint();
  SpatialHelper spatial_helper(labels);

  timer.createTimepoint();
  std::unordered_map<OsmId, PointOfInterest> pois;
  std::transform(labels.begin(),
                 labels.end(),
                 std::inserter(pois, pois.begin()),
                 [](PointOfInterest& poi) {
                   return std::make_pair(poi.get_osm_id(), std::move(poi));
                 });

  // remove dupplicated pois
  for (auto& d : spatial_helper.get_dupplicates()) {
    auto dupplicate = pois.find(d);
    std::cout << "Ignoring dupplicated poi: " << dupplicate->second.print()
              << std::endl;
    pois.erase(dupplicate);
  }

  timer.createTimepoint();
  std::vector<Elimination> result;

  std::priority_queue<Event> Q;

  timer.createTimepoint();
  // initialize
  for (const auto& p : pois) {
    Event evt = predict_collision(p.second, 0., spatial_helper, pois);

    assert(evt.m_evt_type == EventType::UPDATE_EVENT);
    Q.push(evt);
  }

  timer.createTimepoint();
  Time t = 0;
  while (!Q.empty()) {
    auto current_evt = Q.top();
    Q.pop();
    t = current_evt.m_trigger_time;

    auto p1 = pois.find(current_evt.m_coll1);
    if (p1 != pois.end()) {
      if (current_evt.m_evt_type == EventType::UPDATE_EVENT) {
        auto evt = predict_collision(p1->second, t, spatial_helper, pois);
        Q.push(evt);
      } else if (current_evt.m_evt_type == EventType::COLLISION_EVENT) {
        if (current_evt.m_coll1 == current_evt.m_coll2) {
          result.emplace_back(t, p1->second, p1->second);
          break;
        }

        auto p2 = pois.find(current_evt.m_coll2);
        if (p2 != pois.end()) {
          // here p1 and p2 are alive
          if (prefer_p1(p1->second, p2->second, spatial_helper, pois)) {
            //             result.emplace_back(coll_t, std::move(p2->second));
            result.emplace_back(t, p2->second, p1->second);
            spatial_helper.erase(p2->first);
            pois.erase(p2);

            auto evt = predict_collision(p1->second, t, spatial_helper, pois);
            Q.push(evt);
          } else {
            result.emplace_back(t, p1->second, p2->second);
            spatial_helper.erase(p1->first);
            pois.erase(p1);

            auto evt = predict_collision(p2->second, t, spatial_helper, pois);
            Q.push(evt);
          }
        } else {
          // p1 is the only collision partner alive => repredict collision
          auto evt = predict_collision(p1->second, t, spatial_helper, pois);
          Q.push(evt);
        }
      }
    } else {
      // check if p2 is alive and predict its next collision
      // otherwise do nothing
      auto p2 = pois.find(current_evt.m_coll2);
      if (p2 != pois.end()) {
        auto evt = predict_collision(p2->second, t, spatial_helper, pois);
        Q.push(evt);
      };
    }
  }
  timer.stop();

  auto times = timer.getTimes();

  std::cout << "Computation of the elimination order finished." << std::endl;
  std::cout << "Required times in ms were as follows:" << std::endl;
  std::cout << times[0] << "\t Label import" << std::endl;
  std::cout << times[1] << "\t Initialization of the spatial helper"
            << std::endl;
  std::cout << times[2] << "\t Creation of the poi map" << std::endl;
  std::cout << times[4] << "\t Initialization of the algorithm" << std::endl;
  std::cout << times[5] << "\t Main algorithm loop" << std::endl;

  return result;
}

// END class EliminationOrder
}
// END DEFINITION
