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

#include <stb_image.h>
#include <stb_image_write.h>


typedef double PixType;
class Image{
private:
	PixType *data;
	int w, h;
	bool ownsData;
public:

	Image(const std::string &path) : ownsData(true){
		int n;
		// load image as greyscale
		unsigned char *datauc = stbi_load(path.c_str(), &w, &h, &n, 1);

		data = new PixType[w*h];
		for(int i = 0; i < w*h; ++i)
			data[i] = datauc[i];
		
		stbi_image_free(datauc);
	}
	
	int width() const {return w;}
	int height() const {return h;}

	PixType &operator()(int i, int j){
		return data[i*w +j];
	}

	const PixType &operator()(int i, int j) const {
		return data[i*w +j];
	}

	PixType *getData(){
		return data;
	}

	const PixType *getData() const{
		return data;
	}

	void write(const std::string &dest) const {
		std::vector<unsigned char> datauc(w*h);
		for(int i = 0; i < w*h; ++i)
			datauc[i] = (unsigned char)(data[i]);

		stbi_write_png(dest.c_str(), w, h, 1, &datauc[0], w);
	}

	~Image(){
		if(ownsData)
			delete [] data;	
	}

};

struct Cell{
	enum Type {Symbol, Number, List};
    typedef Cell (*proc_type)(const std::vector<Cell> &);
    typedef std::vector<Cell>::const_iterator iter;
    Type type; std::string val; std::vector<Cell> list;
    Cell(Type type = Symbol) : type(type) {}
    Cell(Type type, const std::string & val) : type(type), val(val) {}
};


template <typename EvalReturn> class Visitor{
public:
	typedef std::map<std::string, std::function<EvalReturn (const std::vector<EvalReturn> &)>> 
		FunctionMap;

	typedef std::function<EvalReturn (const std::string &symbol)> SymbolHandler;
	typedef std::function<EvalReturn (const std::string &number)> NumberHandler;

protected:
FunctionMap functionMap;
NumberHandler numberHandler;
SymbolHandler symbolHandler;

public:
	Visitor(){
	}

	EvalReturn eval(const Cell &c){
		switch(c.type){
			case Cell::Number:{
				return numberHandler(c.val.c_str());
			}case Cell::List:{
				std::vector<EvalReturn> evalArgs(c.list.size()-1);
				// eval each argument
				std::transform(c.list.begin()+1, c.list.end(), evalArgs.begin(), 
					[=](const Cell &c) -> EvalReturn{
						return this->eval(c);
					}
				);
				// call function specified by sumbol map with evaled arguments
				return functionMap.at(c.list[0].val)(evalArgs);
			}case Cell::Symbol:{
				if(symbolHandler)
					return symbolHandler(c.val);
				else
					std::runtime_error("Cannot handle symbol.");
			}
			std::runtime_error("Should never get here.");
		}
	}

};

class Calculator : public Visitor<double>{
public:
	Calculator(){
		// standard functions
		functionMap["+"] = [](const std::vector<double> &d){return d[0] + d[1];};
		functionMap["-"] = [](const std::vector<double> &d){return d[0] - d[1];};
		functionMap["/"] = [](const std::vector<double> &d){return d[0] / d[1];};
		functionMap["*"] = [](const std::vector<double> &d){return d[0] * d[1];};
		functionMap["sin"] = [](const std::vector<double> &d){return std::sin(d[0]);};
		functionMap["cos"] = [](const std::vector<double> &d){return std::cos(d[0]);};
		functionMap["tan"] = [](const std::vector<double> &d){return std::tan(d[0]);};
		functionMap["pow"] = [](const std::vector<double> &d){return std::pow(d[0], d[1]);};

		numberHandler = [](const std::string &number){
			return std::atof(number.c_str());
		};


	}

};

class CalculatorFunction : public Calculator{
private:
	std::map<std::string, int> argNameToIndex;
	Cell cell;
public:
	CalculatorFunction(const std::vector<std::string> &names, const Cell &c) : cell(c){
		for(int i = 0; i < names.size(); ++i)
			argNameToIndex[names[i]] = i;
	}

	double operator()(const std::vector<double> &args){
		symbolHandler = [&](const std::string &name) -> double{
			return args[this->argNameToIndex[name]];	
		};

		return eval(cell);
	}
};

	
class CodeGenCalculatorFunction : public Visitor<AsmJit::XmmVar>{
private:
	AsmJit::X86Compiler compiler;
	std::map<std::string, int> argNameToIndex;

	typedef double (*FuncPtrType)(const double * args);
	FuncPtrType generatedFunction;
public:
	CodeGenCalculatorFunction(const std::vector<std::string> &names, const Cell &cell){
		using namespace AsmJit;
		
		functionMap["+"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			// XmmVar resultVar(compiler.newXmmVar());
			compiler.addsd(args[0], args[1]);
			// compiler.movq(resultVar, args[0]);
			return args[0];
		};
	
		functionMap["-"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.subsd(args[0], args[1]);
			return args[0];
		};
	
		functionMap["*"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.mulsd(args[0], args[1]);
			return args[0];
		};
	
		functionMap["/"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.divsd(args[0], args[1]);
			return args[0];
		};


		// functionMap["-"] = [](const std::vector<double> &d){return d[0] - d[1];};
		// functionMap["/"] = [](const std::vector<double> &d){return d[0] / d[1];};
		// functionMap["*"] = [](const std::vector<double> &d){return d[0] * d[1];};

		numberHandler = [&](const std::string &number) -> XmmVar{
			double x = std::atof(number.c_str());
			XmmVar xVar(compiler.newXmmVar());
			SetXmmVar(compiler, xVar, x);
			return xVar;
		};

		for(int i = 0; i < names.size(); ++i)
			argNameToIndex[names[i]] = i;

		symbolHandler = [&](const std::string name) -> XmmVar{
			// Lookup name in args and return AsmJit variable
			// with the arg loaded in.
			// TODO: this could be more efficient - could
			// create one list of XmmVars and use that.
			GpVar ptr(compiler.getGpArg(0));
			XmmVar v(compiler.newXmmVar());
			int offset = argNameToIndex.at(name)*sizeof(double);
			compiler.movsd(v, Mem(ptr, offset));
			return v;
		};

		generatedFunction = generate(cell);
	}

	FuncPtrType generate(const Cell &c){
		compiler.newFunc(AsmJit::kX86FuncConvDefault, 
		                 AsmJit::FuncBuilder1<double, const double *>());
		AsmJit::XmmVar retVar = eval(c);
		compiler.ret(retVar);
		compiler.endFunc();
		return reinterpret_cast<FuncPtrType>(compiler.make());
		
	}

	double operator()(const std::vector<double> &args) const {
		return generatedFunction(&args[0]); 
	}

	private:
	void SetXmmVar(AsmJit::X86Compiler &c, AsmJit::XmmVar &v, double d){
		using namespace AsmJit;
        GpVar tmp(c.newGpVar());
        c.mov(tmp, d);
        c.movq(v, tmp);
        c.unuse(tmp);
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
	/*
	Calculator calc;
    for (;;) {
		CodeGenCalculator cgcalc;
        std::cout << prompt;
        std::string line; std::getline(std::cin, line);
       	Cell cell = read(line);
        std::cout << calc.eval(read(line)) << '\n';
        // std::function<double (void)> f = cgcalc.generate(read(line));
        std::cout << "code gen " << f() << '\n';
    }
    */
}

int main ()
{
	Image im("/home/luke/projects/pixslam/example_data/lena.png");
	im.write("test.png");
	
	// repl(">");
	std::vector<std::string> argNames = {"x", "y", "z"};
	std::vector<double> args = {1.5, 2.5};
	std::string functionCode = ("(+ x (/ x (- y x)))");
	Cell functionCell = read(functionCode);
	CalculatorFunction function(argNames, functionCell);
	CodeGenCalculatorFunction cgFunction(argNames, functionCell);
	std::cout << "Interpreted output: " << function(args) << std::endl;
	std::cout << "Code gen output: " << cgFunction(args) << std::endl;

	int repetitions = 1000000;
	for(int i = 0; i < repetitions; ++i)
		cgFunction(args);
	
	return 0;
}

