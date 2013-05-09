#pragma once

#include <vector>
#include <string>

namespace pixslam{

struct Cell{
    enum Type {Symbol, Number, List};
    typedef Cell (*proc_type)(const std::vector<Cell> &);
    typedef std::vector<Cell>::const_iterator iter;
    Type type; std::string val; std::vector<Cell> list;
    Cell(Type type = Symbol) : type(type) {}
    Cell(Type type, const std::string & val) : type(type), val(val) {}
};

// convert given Cell to a Lisp-readable string
std::string to_string(const Cell & exp);

// return the Lisp expression represented by the given string
Cell read(const std::string & s);

}
