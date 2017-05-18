#include <iostream>
#include <vector>

#include "datastorage.h"

int main ( int argc, char** argv )
{
  using DataStorage = growing_balls::DataStorage;
  using Element = growing_balls::DataStorage::Element;
  
    std::vector<Element> input_elems;
    input_elems.emplace_back(0., 0., 0);
    input_elems.emplace_back(0., 1., 1);
    input_elems.emplace_back(.5, 1., 2);
    input_elems.emplace_back(.5, -2., 3);
    
    std::cout << "Bulk inserting " << input_elems.size() << " elements ..." << std::flush;
    DataStorage storage;
    storage.insert(input_elems.begin(), input_elems.end());
    std::cout << "\t\tsuccess!" << std::endl;
    
    std::cout << "Inserting singĺe element " << std::flush;
    Element e = Element(0., 1., 4);
    auto id = storage.insert(e);
    std::cout << "\t\tsuccessfully inserted element with id " << id << std::endl;
    
    std::cout << "Requesting adjacents of element with id " << id << " ..." << std::endl; 
    auto visit =  [](const Element& elem) {
      std::cout << "Element: " << elem.get_info() << " got id " << elem.get_id() << std::endl;
    };
    storage.neighborhood(id, visit);
    
    std::cout << "\tfinished." << std::endl;
}
