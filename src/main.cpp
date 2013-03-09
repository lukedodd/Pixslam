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
			datauc[i] = (unsigned char)(std::max(std::min(data[i], 255.0), 0.0));

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

				return functionHandler(c.list[0].val, evalArgs);
			}case Cell::Symbol:{
				return symbolHandler(c.val);
			}
			throw std::runtime_error("Should never get here.");
		}
	}

	virtual ~Visitor(){}

protected:
	virtual EvalReturn symbolHandler(const std::string &symbol) = 0;
	virtual EvalReturn functionHandler(const std::string &functionName, 
	                                   const std::vector<EvalReturn> &args) = 0;
	virtual EvalReturn numberHandler(const std::string &number) = 0;

};

class Calculator : public Visitor<double>{
	typedef std::function<double (const std::vector<double> &)> BuiltInFunctionHandler;
	std::map<std::string, BuiltInFunctionHandler> functionHandlerMap;
public:
	Calculator(){
		// standard functions
		functionHandlerMap["+"] = [](const std::vector<double> &d){return d[0] + d[1];};
		functionHandlerMap["-"] = [](const std::vector<double> &d){return d[0] - d[1];};
		functionHandlerMap["/"] = [](const std::vector<double> &d){return d[0] / d[1];};
		functionHandlerMap["*"] = [](const std::vector<double> &d){return d[0] * d[1];};
		functionHandlerMap["sin"] = [](const std::vector<double> &d){return std::sin(d[0]);};
		functionHandlerMap["cos"] = [](const std::vector<double> &d){return std::cos(d[0]);};
		functionHandlerMap["tan"] = [](const std::vector<double> &d){return std::tan(d[0]);};
		functionHandlerMap["pow"] = [](const std::vector<double> &d){return std::pow(d[0], d[1]);};

	}

	virtual ~Calculator(){}

protected:

	virtual double functionHandler(const std::string &functionName, 
	                               const std::vector<double> &args){
		return functionHandlerMap.at(functionName)(args);
	}

	virtual double numberHandler(const std::string &number){
		return std::atof(number.c_str());
	}

	virtual double symbolHandler(const std::string &symbol){
		throw std::runtime_error("Cannot handle symbol.");
	}

};

class CalculatorFunction : public Calculator{
private:
	std::map<std::string, int> argNameToIndex;
	std::vector<double> inputArgs;
	Cell cell;
public:
	CalculatorFunction(const std::vector<std::string> &names, const Cell &c) : cell(c){
		for(int i = 0; i < names.size(); ++i)
			argNameToIndex[names[i]] = i;
	}

	double operator()(const std::vector<double> &args){
		inputArgs = args;
		return eval(cell);
	}

	virtual ~CalculatorFunction(){}

protected:
	virtual double symbolHandler(const std::string &symbol){
		return inputArgs[this->argNameToIndex[symbol]];	
	}
};

	
class CodeGenCalculatorFunction : public Visitor<AsmJit::XmmVar>{
private:
	typedef std::function<AsmJit::XmmVar (const std::vector<AsmJit::XmmVar> &)> 
	        BuiltInFunctionHandler;
	std::map<std::string, BuiltInFunctionHandler> functionHandlerMap;

	AsmJit::X86Compiler compiler;
	std::map<std::string, int> argNameToIndex;

	typedef double (*FuncPtrType)(const double * args);
	FuncPtrType generatedFunction;
public:
	CodeGenCalculatorFunction(const std::vector<std::string> &names, const Cell &cell){
		using namespace AsmJit;
		
		functionHandlerMap["+"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.addsd(args[0], args[1]);
			return args[0];
		};
	
		functionHandlerMap["-"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.subsd(args[0], args[1]);
			return args[0];
		};
	
		functionHandlerMap["*"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.mulsd(args[0], args[1]);
			return args[0];
		};
	
		functionHandlerMap["/"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
			compiler.divsd(args[0], args[1]);
			return args[0];
		};


		for(int i = 0; i < names.size(); ++i)
			argNameToIndex[names[i]] = i;

	
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

	virtual ~CodeGenCalculatorFunction(){}

protected:

	virtual AsmJit::XmmVar functionHandler(const std::string &functionName, 
	                                       const std::vector<AsmJit::XmmVar> &args){
		return functionHandlerMap.at(functionName)(args);
	}

	virtual AsmJit::XmmVar numberHandler(const std::string &number){
		double x = std::atof(number.c_str());
		AsmJit::XmmVar xVar(compiler.newXmmVar());
		SetXmmVar(compiler, xVar, x);
		return xVar;
	};

	virtual AsmJit::XmmVar symbolHandler(const std::string &name){
			// Lookup name in args and return AsmJit variable
			// with the arg loaded in.
			// TODO: this could be more efficient - could
			// create one list of XmmVars and use that.
			using namespace AsmJit;
			GpVar ptr(compiler.getGpArg(0));
			XmmVar v(compiler.newXmmVar());
			int offset = argNameToIndex.at(name)*sizeof(double);
			compiler.movsd(v, Mem(ptr, offset));
			return v;
	};

private:
	void SetXmmVar(AsmJit::X86Compiler &c, AsmJit::XmmVar &v, double d){
		using namespace AsmJit;
		// I thought this was better, but it didn't actually work...
        // GpVar tmp(c.newGpVar());
        // c.mov(tmp, d);
        // c.movq(v, tmp);
        // c.unuse(tmp);

		GpVar half(c.newGpVar());
        uint64_t *i = reinterpret_cast<uint64_t*>(&d);
        c.mov(half, i[0]);
        c.movq(v, half);
        c.unuse(half);
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

int main (int argc, char *argv[])
{
	std::string outputImage = "out.png";

	if(argc < 2){
		std::cout << "Please supply an input image! \n";
		return 1;
	}
	
	if(argc >= 3)
		outputImage = argv[2];

	Image im(argv[1]);
	
	// repl(">");
	std::vector<std::string> argNames(1, "x");
	std::vector<double> args(1, 1.5);
	std::string functionCode = "(* x 2)";
	Cell functionCell = read(functionCode);

	CalculatorFunction function(argNames, functionCell);
	CodeGenCalculatorFunction cgFunction(argNames, functionCell);

	std::cout << "Interpreted output: " << function(args) << std::endl;
	std::cout << "Code gen output: " << cgFunction(args) << std::endl;

	for(int i = 0; i < im.width(); ++i){
		for(int j = 0; j < im.height(); ++j){
			args[0] = im(i,j);
			im(i,j) = cgFunction(args);
			// im(i,j) = function(args);
		}
	}

	im.write(outputImage);
	
	return 0;
}

