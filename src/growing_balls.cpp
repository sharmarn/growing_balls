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
    auto et = e.first;
    auto poi = e.second;
    std::cout << et << "\t" << poi.print() << std::endl;
  }
}
