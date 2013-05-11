#include "Parser.h"

#include <list>

namespace pixslam{

namespace{
    // Anonymous nested name space for functions required to implement parsing/printing
    bool isWhiteSpace(char c){
        return c == ' ' || c == '\t' || c == '\n';
    }

    bool isDig(char c) { return isdigit(static_cast<unsigned char>(c)) != 0; }

    // convert given string to list of tokens
    // originally from: 
    // http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
    std::list<std::string> tokenize(const std::string & str){
        std::list<std::string> tokens;
        const char * s = str.c_str();
        while (*s) {
            while(isWhiteSpace(*s)) // ignore whitespace
                ++s;

            if(*s == ';'){ 
                // skip to newline after a comment
                while(*s != '\n')
                    ++s; 
            }else if(*s == '(' || *s == ')'){
                tokens.push_back(*s++ == '(' ? "(" : ")");
            }else{
                const char * t = s;
                while(*t && !isWhiteSpace(*t) && *t != '(' && *t != ')')
                    ++t;
                tokens.push_back(std::string(s, t));
                s = t;
            }
        }
        return tokens;
    }

    // numbers become Numbers; every other token is a Symbol
    // originally from: 
    // http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
    Cell atom(const std::string & token)
    {
        if (isDig(token[0]) || (token[0] == '-' && isDig(token[1])))
            return Cell(Cell::Number, token);
        return Cell(Cell::Symbol, token);
    }
    // return the Lisp expression in the given tokens
    // originally from: 
    // http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
    Cell readFrom(std::list<std::string> & tokens)
    {
        const std::string token(tokens.front());
        tokens.pop_front();
        if (token == "(") {
            Cell c(Cell::List);
            while (tokens.front() != ")")
                c.list.push_back(readFrom(tokens));
            tokens.pop_front();
            return c;
        }
        else
            return atom(token);
    }
}

// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
Cell cellFromString(const std::string & s)
{
    std::list<std::string> tokens(tokenize(s));
    return readFrom(tokens);
}

// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
std::string cellToString(const Cell & exp)
{
    if (exp.type == Cell::List) {
        std::string s("(");
        for (Cell::iter e = exp.list.begin(); e != exp.list.end(); ++e)
            s += cellToString(*e) + ' ';
        if (s[s.size() - 1] == ' ')
            s.erase(s.size() - 1);
        return s + ')';
    }
    return exp.val;
}

}

