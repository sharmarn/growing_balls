#include <iostream>
#include <vector>

#include "datastorage.h"
#include "eliminationorder.h"
#include "spatialhelper.h"
#include "textinput.h"

namespace {
  using TextInput = growing_balls::TextInput;
  using OsmId = TextInput::OsmId;
  using PointOfInterest = growing_balls::TextInput::PointOfInterest;
  using DataStorage = growing_balls::DataStorage<OsmId>;
  using LabelElement = DataStorage::Element;

LabelElement
element_from(PointOfInterest poi)
{
  return LabelElement(poi.get_lat(), poi.get_lon(), poi.get_osm_id());
}
}

int
main(int argc, char** argv)
{ 
  if (argc < 2) {
    std::cout
      << "Please provide a valid .complete.txt input file and a valid 'osm_id'!"
      << std::endl;
    return 1;
  }
  
//   std::cout << "Inserting data to the spatial helper ..." << std::endl;
//   growing_balls::SpatialHelper sph(TextInput::import_label(argv[1]));
//   OsmId query = 10;
//   std::cout << "NearestNeighbor of element #" << query << " ..." << std::endl;
//   std::cout << "\thas id #" << sph.get_nearest_neighbor(query) << std::endl;
//   growing_balls::SpatialHelper::Distance d = 22300000;
//   std::cout << "Requesting elements in " << d << "-neighborhood ..." << std::endl;
//   for (auto& e : sph.get_in_range(query, d)) {
//     std::cout << e << std::endl;
//   }
//   
//   std::cout << "Inserting data to the data storage ..." << std::endl;
//   DataStorage ds;
//   OsmId id;
//   
//   for (auto& poi : TextInput::import_label(argv[1])) {
//     id = ds.insert(element_from(poi));
//   }
//   std::cout << "Outputting the neighbor elements of i_id ..." << id << std::endl;
//   auto visitor = [](const LabelElement& elem) {
//     std::cout << "Element #" << elem.get_info() << " at (" << elem.get_coord_1() << ", " << elem.get_coord_2() << ") with internal id " << elem.get_id() << std::endl;
//   };
//   ds.visit_neighborhood(id, visitor);
//  
//   return 0;
  
  growing_balls::EliminationOrder eo(argv[1]);
  auto es = eo.compute_elimination_order();
  
  for (auto& e : es) {
    auto et = e.first;
    auto poi = e.second;
    std::cout << et << "\t" << poi.print() << std::endl;
  }
}
