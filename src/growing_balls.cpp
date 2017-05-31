#include <iostream>
#include <string>
#include <vector>

#include "eliminationorder.h"
#include "io.h"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Please provide a valid .complete.txt input file!"
              << std::endl;
    return 1;
  }

  growing_balls::EliminationOrder eo;
  auto es = eo.compute_elimination_order(argv[1]);

  std::vector<growing_balls::IO::PointOfInterest> elimination_order;
  for (auto &e : es) {
    auto et = e.m_elimination_time;
    auto poi = e.m_eliminated;
    auto elim_p = e.m_eliminated_by;

    poi.set_elimination(et, elim_p.get_osm_id());
    elimination_order.push_back(std::move(poi));
  }

  std::string out_path = "elim_out.raw";
  growing_balls::IO::export_eliminationorder(out_path, elimination_order);
}
