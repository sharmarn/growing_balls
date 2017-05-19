#include <iostream>
#include <vector>

#include "datastorage.h"
#include "spatialfunctions.h"
#include "textinput.h"

namespace {
using LabelElement =
  growing_balls::DataStorage<growing_balls::TextInput::OsmId>::Element;

LabelElement
element_from(growing_balls::TextInput::PointOfInterest poi)
{
  return LabelElement(poi.get_lat(), poi.get_lon(), poi.get_osm_id());
}
}

int
main(int argc, char** argv)
{
  growing_balls::TextInput::PointOfInterest p("48.00000000000000000 "
                                              "11.90000000000000036 0 180.000 "
                                              "25873009 'Bayern' 1.000");

  std::cout << p.print() << std::endl;

  auto labels = growing_balls::TextInput::import_label("test.complete.txt");

  std::cout << "Imported " << labels.size() << " elements" << std::endl;

  std::vector<LabelElement> lbl_elems;
  for (auto& lbl : labels) {
    lbl_elems.push_back(element_from(lbl));
  }

  growing_balls::DataStorage<growing_balls::TextInput::OsmId> store;
  store.insert(lbl_elems.begin(), lbl_elems.end());

  return 0;

  using DataStorage = growing_balls::DataStorage<std::size_t>;
  using Element = growing_balls::DataStorage<std::size_t>::Element;

  std::vector<Element> input_elems;
  input_elems.emplace_back(0., 0., 0);
  input_elems.emplace_back(0., 1., 1);
  input_elems.emplace_back(.5, 1., 2);
  input_elems.emplace_back(.5, -2., 3);

  std::cout << "Bulk inserting " << input_elems.size() << " elements ..."
            << std::flush;
  DataStorage storage;
  storage.insert(input_elems.begin(), input_elems.end());
  std::cout << "\t\tsuccess!" << std::endl;

  std::cout << "Inserting singĺe element " << std::flush;
  Element e = Element(0., 1., 4);
  auto id = storage.insert(e);
  std::cout << "\t\tsuccessfully inserted element with id " << id << std::endl;

  std::cout << "Requesting adjacents of element with id " << id << " ..."
            << std::endl;
  auto visit = [](const Element& elem) {
    std::cout << "Element: " << elem.get_info() << " got id " << elem.get_id()
              << std::endl;
  };
  storage.neighborhood(id, visit);

  std::cout << "\tfinished." << std::endl;
}
