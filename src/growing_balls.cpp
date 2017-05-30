#include <iostream>
#include <vector>

#include "eliminationorder.h"

int
main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Please provide a valid .complete.txt input file!"
              << std::endl;
    return 1;
  }

  growing_balls::EliminationOrder eo(argv[1]);
  auto es = eo.compute_elimination_order(argv[1]);

  for (auto& e : es) {
    auto et = e.m_elimination_time;
    auto poi = e.m_eliminated;
    auto elim_p = e.m_eliminated_by;
    std::cout << et << "\t\t" << poi.print() << "\t\t" << elim_p.get_priority()
              << std::endl;
  }
}
