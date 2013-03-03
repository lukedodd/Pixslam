#include <vector>
#include <list>
#include <iostream>
#include <map>
#include <functional>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#include <asmjit/asmjit.h>



struct Cell{
	enum Type {Symbol, Number, List};
    typedef Cell (*proc_type)(const std::vector<Cell> &);
    typedef std::vector<Cell>::const_iterator iter;
    Type type; std::string val; std::vector<Cell> list;
    Cell(Type type = Symbol) : type(type) {}
    Cell(Type type, const std::string & val) : type(type), val(val) {}
};


class Visitor{
private:
typedef std::map<std::string, std::function<double (const std::vector<double> &)>> SymbolMap;
SymbolMap symbolMap;
public:
	Visitor(){
		symbolMap["+"] = [](const std::vector<double> &d){return d[0] + d[1];};
		symbolMap["-"] = [](const std::vector<double> &d){return d[0] - d[1];};
		symbolMap["/"] = [](const std::vector<double> &d){return d[0] / d[1];};
		symbolMap["*"] = [](const std::vector<double> &d){return d[0] * d[1];};
		symbolMap["sin"] = [](const std::vector<double> &d){return std::sin(d[0]);};
		symbolMap["cos"] = [](const std::vector<double> &d){return std::cos(d[0]);};
		symbolMap["tan"] = [](const std::vector<double> &d){return std::tan(d[0]);};
		symbolMap["pow"] = [](const std::vector<double> &d){return std::pow(d[0], d[1]);};
	}

	double eval(const Cell &c){
		switch(c.type){
			case Cell::Number:{
				return std::atof(c.val.c_str());
				break;
			}case Cell::List:{
				std::vector<double> evalArgs(c.list.size()-1);
				std::transform(c.list.begin()+1, c.list.end(), evalArgs.begin(), 
					[=](const Cell &c) -> double{
						return this->eval(c);
					}
				);
				return symbolMap.at(c.list[0].val)(evalArgs);
			}case Cell::Symbol:{
				std::runtime_error("Symbol not expected.");
			}
			return 0.0;
		}
	}
};


// convert given string to list of tokens
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
std::list<std::string> tokenize(const std::string & str){
    std::list<std::string> tokens;
    const char * s = str.c_str();
    while (*s) {
        while (*s == ' ')
            ++s;
        if (*s == '(' || *s == ')')
            tokens.push_back(*s++ == '(' ? "(" : ")");
        else {
            const char * t = s;
            while (*t && *t != ' ' && *t != '(' && *t != ')')
                ++t;
            tokens.push_back(std::string(s, t));
            s = t;
        }
    }
    return tokens;
}

bool isdig(char c) { return isdigit(static_cast<unsigned char>(c)) != 0; }

// numbers become Numbers; every other token is a Symbol
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
Cell atom(const std::string & token)
{
    if (isdig(token[0]) || (token[0] == '-' && isdig(token[1])))
        return Cell(Cell::Number, token);
    return Cell(Cell::Symbol, token);
}

// return the Lisp expression in the given tokens
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
Cell read_from(std::list<std::string> & tokens)
{
    const std::string token(tokens.front());
    tokens.pop_front();
    if (token == "(") {
        Cell c(Cell::List);
        while (tokens.front() != ")")
            c.list.push_back(read_from(tokens));
        tokens.pop_front();
        return c;
    }
    else
        return atom(token);
}

// return the Lisp expression represented by the given string
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
Cell read(const std::string & s)
{
    std::list<std::string> tokens(tokenize(s));
    return read_from(tokens);
}

// convert given Cell to a Lisp-readable string
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
std::string to_string(const Cell & exp)
{
    if (exp.type == Cell::List) {
        std::string s("(");
        for (Cell::iter e = exp.list.begin(); e != exp.list.end(); ++e)
            s += to_string(*e) + ' ';
        if (s[s.size() - 1] == ' ')
            s.erase(s.size() - 1);
        return s + ')';
    }
    return exp.val;
}

// read-eval-print-loop
void repl(const std::string & prompt )
{
	Visitor v;
    for (;;) {
        std::cout << prompt;
        std::string line; std::getline(std::cin, line);
        std::cout << to_string(read(line)) << '\n';
        std::cout << v.eval(read(line)) << '\n';
    }
}

int main ()
{
	repl(">");
}

