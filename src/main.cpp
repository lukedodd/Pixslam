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
#include <numeric>
#include <sstream>
#include <fstream>

#include <asmjit/asmjit.h>

#include <stb_image.h>
#include <stb_image_write.h>


typedef double PixType;
class Image{
private:
    PixType *data;
    int w, h;
    int s;
    bool ownsData;
public:

    Image(const std::string &path) : ownsData(true){
        int n;
        // load image as greyscale
        unsigned char *datauc = stbi_load(path.c_str(), &w, &h, &n, 1);
        s = w;

        data = new PixType[w*h];
        for(int i = 0; i < w*h; ++i)
            data[i] = datauc[i]/255.0;

        stbi_image_free(datauc);
    }

    Image(int w, int h) : w(w), h(h), s(w), ownsData(true){
        data = new PixType[w*h];
        std::fill(data, data+(w*h), 0.0);
    }

    Image(int w_, int h_, int s_) : w(w_), h(h_), s(s_), ownsData(true){
        data = new PixType[h*s];
        std::fill(data, data+(h*s), 0.0);
    }



    Image(PixType *data, int w, int h, int s)
        : data(data), w(w), h(h), s(s), ownsData(false)
    {
    }

    Image(const Image &original, int padx, int pady)
        : w(original.width() + padx*2), h(original.height() + pady*2), s(w), ownsData(true){
        data = new PixType[w*h];
        std::fill(data, data+(w*h), 0.0);
        for(int i = 0; i < original.height(); ++i){
            for(int j = 0; j < original.width(); ++j){
                (*this)(i+pady, j+padx) = original(i,j);
            }
        }
    }

    // Forbid copy and assignment for now, allow move.
    Image(const Image &) = delete;
    Image &operator=(const Image&) = delete;

    Image(Image&& other) : 
        data(other.data), w(other.w), h(other.h), s(other.s), ownsData(other.ownsData){
        other.data = 0;
    }

    int width() const {return w;}
    int height() const {return h;}
    int stride() const {return s;}

    PixType &operator()(int i, int j){
        return data[i*s +j];
    }

    const PixType &operator()(int i, int j) const {
        return data[i*s +j];
    }

    PixType *getData(){
        return data;
    }

    const PixType *getData() const{
        return data;
    }

    void write(const std::string &dest) const {
        std::vector<unsigned char> datauc(w*h);

        for(int i = 0; i < w; ++i)
            for(int j = 0; j < h; ++j)
                datauc[i*w + j] = (unsigned char)(std::max(
                                  std::min((*this)(i,j)*255.0, 255.0), 0.0));

        stbi_write_png(dest.c_str(), w, h, 1, &datauc[0], w);
    }

    ~Image(){
        if(ownsData && data)
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
        }
        throw std::runtime_error("Should never get here.");
    }

    virtual ~Visitor(){}

protected:
    virtual EvalReturn symbolHandler(const std::string &symbol) = 0;
    virtual EvalReturn functionHandler(const std::string &functionName, 
            const std::vector<EvalReturn> &args) = 0;
    virtual EvalReturn numberHandler(const std::string &number) = 0;

};

class CodeGenCalculatorFunction : public Visitor<AsmJit::XmmVar>{
private:
    typedef std::function<AsmJit::XmmVar (const std::vector<AsmJit::XmmVar> &)> 
        BuiltInFunctionHandler;

    typedef const double * const * Arguments;
    std::map<std::string, BuiltInFunctionHandler> functionHandlerMap;

    AsmJit::X86Compiler compiler;
    std::map<std::string, size_t> argNameToIndex;

    typedef void (*FuncPtrType)(Arguments, size_t, size_t, size_t, double *);
    FuncPtrType generatedFunction;

    AsmJit::GpVar currentI;
    AsmJit::GpVar currentJ;
    AsmJit::GpVar currentIndex;

    AsmJit::GpVar w;
    AsmJit::GpVar h;
    AsmJit::GpVar stride;
    AsmJit::GpVar n;
    AsmJit::GpVar out;

    AsmJit::XmmVar wd;
    AsmJit::XmmVar hd;

    AsmJit::XmmVar zero;
    AsmJit::XmmVar one;
    std::vector<AsmJit::GpVar> argv;

    std::map<std::string, AsmJit::XmmVar> symbols;

    AsmJit::FileLogger logger;

public:
    CodeGenCalculatorFunction(const Cell &cell) : logger(stdout) {

        // compiler.setLogger(&logger);
       
        // Check cell is of form ((arg list) (expr))
        if(!(cell.type == Cell::List && cell.list.size() == 2 &&
             cell.list[0].type == Cell::List && cell.list[1].type == Cell::List))
           throw std::runtime_error("Function cell must be of form ((arg1 arg2 ...) (code))");


        const Cell &argsCell = cell.list[0];
        const Cell &code = cell.list[1];

        // Load arguments
        std::vector<std::string> argNames;
        for(Cell c : argsCell.list){
            if(c.type == Cell::Symbol)
                argNames.push_back(c.val);
            else
                throw std::runtime_error(
                    "Function cell must be of form ((arg1 arg2 ...) (code))");
        }
        

        for(size_t i = 0; i < argNames.size(); ++i)
            argNameToIndex[argNames[i]] = i;


        // Specify how to deal with function calls.
        PopulateBuiltInFunctionHandlerMap();
        // Generate the code!
        generatedFunction = generate(code);
    }

    void operator()(const std::vector<Image> &images, Image &out) const {
        if(images.empty())
            throw std::runtime_error("must have at least one input image.");

        // check all images have the same dimension + stride.
        for(size_t i = 1; i < images.size(); ++i)
            if(   images[i-1].width() != images[i].width()
               || images[i-1].height() != images[i].height()
               || images[i-1].stride() != images[i].stride())
                throw std::runtime_error("all input images must have same dimensions.");

        if(   images[0].width() != out.width()
           || images[0].height() != out.height()
           || images[0].stride() != out.stride())
            throw std::runtime_error("all input images and output must have same dimensions.");


        std::vector<const double*> dataPtrs;
        for(const Image &im : images)
            dataPtrs.push_back(im.getData());

        generatedFunction(&dataPtrs[0], images[0].width(), images[0].height(), images[0].stride(),
                          out.getData()); 
    }

    // "Lower level" call for image data from other sources (e.g. opencv)
    void operator()(const std::vector<const double *> &args, 
                    int w, int h, int stride,
                    double *out) const {
        generatedFunction(&args[0], w, h, stride, out); 
    }


    virtual size_t getNumArgs() const {return argNameToIndex.size();}

    virtual ~CodeGenCalculatorFunction(){
        AsmJit::MemoryManager::getGlobal()->free((void*)generatedFunction); 
    }


private:
    FuncPtrType generate(const Cell &c){
        using namespace AsmJit;
        compiler.newFunc(AsmJit::kX86FuncConvDefault, 
                AsmJit::FuncBuilder5<void, Arguments, size_t, size_t, size_t, double *>());


        GpVar pargv = compiler.getGpArg(0);
        for(size_t i = 0; i < argNameToIndex.size(); ++i){
            argv.push_back(compiler.newGpVar());
            compiler.mov(argv.back(), ptr(pargv, i*sizeof(double)));
        }

        zero = compiler.newXmmVar();
        one = compiler.newXmmVar();
        SetXmmVar(compiler, zero, 0.0);
        SetXmmVar(compiler, one, 1.0);

        w = compiler.getGpArg(1);
        h = compiler.getGpArg(2);
        stride = compiler.getGpArg(3);
        out = compiler.getGpArg(4);

        wd = compiler.newXmmVar();
        hd = compiler.newXmmVar();
        compiler.cvtsi2sd(wd, w);
        compiler.cvtsi2sd(hd, h);
        symbols["w"] = wd;
        symbols["h"] = hd;


        // Perpare loop vars
        n = compiler.newGpVar();
        compiler.mov(n, w);
        compiler.imul(n, h);
        currentIndex = compiler.newGpVar();
        compiler.mov(currentIndex, imm(0));

        currentI = compiler.newGpVar();
        currentJ = compiler.newGpVar();
        compiler.mov(currentI, imm(0));
        compiler.mov(currentJ, imm(0));

        // for i = 0..h
        // for j = 0..w
        Label startLoop(compiler.newLabel());
        compiler.bind(startLoop);
        {
            compiler.mov(currentIndex, currentI);
            compiler.imul(currentIndex, stride);
            compiler.add(currentIndex, currentJ);
            // im(i,j) = f(x)
            AsmJit::XmmVar retVar = eval(c);
            compiler.movq(ptr(out, currentIndex, kScale8Times), retVar);

        }
        compiler.add(currentJ, imm(1));
        compiler.cmp(currentJ, w);
        compiler.jne(startLoop);
        compiler.mov(currentJ, imm(0));
        compiler.add(currentI, imm(1));
        compiler.cmp(currentI, h);
        compiler.jne(startLoop);


        compiler.endFunc();
        return reinterpret_cast<FuncPtrType>(compiler.make());
    }

private:

    virtual AsmJit::XmmVar functionHandler(const std::string &functionName, 
            const std::vector<AsmJit::XmmVar> &args){
        using namespace AsmJit;

        // try builtin function lookup first
        auto it = functionHandlerMap.find(functionName);
        if(it != functionHandlerMap.end()){
            return it->second(args);
        }


        if(functionName[0] == '@'){
            std::string imageName = std::string(functionName.begin()+1, functionName.end());

            if(argNameToIndex.find(imageName) == argNameToIndex.end())
                std::runtime_error("Absolute indexing with unknown image " + imageName);

            // Absolute indexing
            GpVar i = compiler.newGpVar();
            GpVar j = compiler.newGpVar();
            compiler.cvtsd2si(i, args[0]);
            compiler.cvtsd2si(j, args[1]);

            GpVar index = compiler.newGpVar();
            compiler.mov(index, i);
            compiler.imul(index, stride);
            compiler.add(index, j);

            GpVar pImage = argv[argNameToIndex.at(imageName)];
            XmmVar v(compiler.newXmmVar());
            compiler.movsd(v, ptr(pImage, index, kScale8Times));
            return v;
        }else{
            std::string imageName = functionName;

            if(argNameToIndex.find(imageName) == argNameToIndex.end())
                std::runtime_error("Unknown function (or image) " + imageName);

            // Otherwise they must have been doing an image lookup...
            GpVar i = compiler.newGpVar();
            GpVar j = compiler.newGpVar();
            compiler.mov(i, currentI);
            compiler.mov(j, currentJ);

            GpVar iOffset = compiler.newGpVar();
            GpVar jOffset = compiler.newGpVar();
            
            // Convert double to int for indexing
            // TODO: figure out a better way.
            compiler.cvtsd2si(iOffset, args[0]);
            compiler.cvtsd2si(jOffset, args[1]);

            compiler.add(i, iOffset);
            compiler.add(j, jOffset);

            GpVar index = compiler.newGpVar();
            compiler.mov(index, i);
            compiler.imul(index, stride);
            compiler.add(index, j);

            
            GpVar pImage = argv[argNameToIndex.at(imageName)];
            XmmVar v(compiler.newXmmVar());
            compiler.movsd(v, ptr(pImage, index, kScale8Times));
            return v;
        }
               
    }

    virtual AsmJit::XmmVar numberHandler(const std::string &number){
        double x = std::atof(number.c_str());
        AsmJit::XmmVar xVar(compiler.newXmmVar());
        SetXmmVar(compiler, xVar, x);
        return xVar;
    };

    virtual AsmJit::XmmVar symbolHandler(const std::string &name){
        using namespace AsmJit;
        // Use of an argument as a symbol not a function call 
        // is equivanlent to x_i_j
        if(argNameToIndex.find(name) != argNameToIndex.end()){
            GpVar pImage = argv[argNameToIndex.at(name)];
            XmmVar v(compiler.newXmmVar());
            compiler.movsd(v, ptr(pImage, currentIndex, kScale8Times));
            return v;
        }else if(name == "i" || name == "j"){ // special symbols
            XmmVar v(compiler.newXmmVar());
            AsmJit::GpVar index = name == "i" ? currentI : currentJ;
            compiler.cvtsi2sd(v, index);
            return v;
        }else if(symbols.find(name) != symbols.end()){
            return symbols[name];
        }
        
        throw std::runtime_error("Unable to find symbol: " + name);
    };

    void SetXmmVar(AsmJit::X86Compiler &c, AsmJit::XmmVar &v, double d){
        // TODO: store constants in memory and load from there...
        using namespace AsmJit;
        GpVar gpreg(c.newGpVar());
        uint64_t *i = reinterpret_cast<uint64_t*>(&d);
        c.mov(gpreg, i[0]);
        c.movq(v, gpreg);
        c.unuse(gpreg);
    }

    void PopulateBuiltInFunctionHandlerMap(){
        using namespace AsmJit;


        // TODO: This could be cleaner, the pattern is the same up to max,
        // but we're not sharing any code.
        functionHandlerMap["+"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.addsd(ret, a);
            });
            return ret;
        };


        functionHandlerMap["-"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.subsd(ret, a);
            });
            return ret;
        };


        functionHandlerMap["*"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.mulsd(ret, a);
            });
            return ret;
        };


        functionHandlerMap["/"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.divsd(ret, a);
            });
            return ret;
        };


        functionHandlerMap["min"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.minsd(ret, a);
            });
            return ret;
        };



        functionHandlerMap["max"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);
            std::for_each(args.begin()+1, args.end(),  [&](const XmmVar &a){
                compiler.maxsd(ret, a);
            });
            return ret;
        };


        functionHandlerMap["<"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);

            compiler.cmpsd(ret, args[1], 1);
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap[">"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 1);
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap["<="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);

            compiler.cmpsd(ret, args[1], 2);
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap[">="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 2);
            compiler.andpd(ret, one);
            return ret;
        };
 
    }

};


// convert given string to list of tokens
// originally from: 
// http://howtowriteaprogram.blogspot.co.uk/2010/11/lisp-interpreter-in-90-lines-of-c.html
std::list<std::string> tokenize(const std::string & str){
    std::list<std::string> tokens;
    const char * s = str.c_str();
    while (*s) {
        while(*s == ' ' || *s == '\t' || *s == '\n') // ignore whitespace
            ++s;

        if(*s == ';'){ // skip to newline after a comment
            while(*s != '\n')
                ++s;
        }else if(*s == '(' || *s == ')'){
            tokens.push_back(*s++ == '(' ? "(" : ")");
        }else{
            const char * t = s;
            while(*t && *t != ' ' && *t != '(' && *t != ')')
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

int main (int argc, char *argv[])
{

    if(argc < 3){
        std::cout << "Usage:\n\n";
        std::cout << "    pixslam <code> <input-images> <output>\n\n";
        std::cout << "Code can either be supplied directly, or as a file path to read in.\n";
        std::cout << "The number of input images read is depndent on the supplied code.\n";
        std::cout << "The output argument is optional, defaults to out.png.\n\n";
        std::cout << "e.g:\n";
        std::cout << "Muliply image by 2 and output to out.png.\n\n";
        std::cout << "    pixslam \"((A) (* A 2))\" image.png\n\n";
        std::cout << "If file mult_by_two.pixslam contains \"(* A 2)\" then the following \n";
        std::cout << "multiplies image.png by 2 and output to image_times_two.png.\n\n";
        std::cout << "    pixslam mult_by_two.pixslam image.png image_times_two.png\n\n";
        std::cout << "Blend two images together equally and output to blend.png.\n";
        std::cout << "    pixslam ((A B) (* 0.5 (+ A B))) image1.png image2.png blend.png\n\n";
        return 1;
    }

    // See if first arg is a file and read code from it.
    std::string codeString;
    std::ifstream ifs(argv[1]);
    if(ifs){
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        codeString = buffer.str();
    }

    if(codeString.empty()) // If that didn't work interpret first arg as code.
        codeString = argv[1];

    // Generate code.
    Cell code = read(codeString);
    CodeGenCalculatorFunction cgFunction(code);

    // Read image from second arg
    int padding = 5;
    std::vector<Image> inputImages;
    for(size_t i = 0; i < cgFunction.getNumArgs(); ++i){
        Image im(argv[2+i]);
        inputImages.emplace_back(im, padding, padding);
    }

    // Remaining arg, if preset is our output destination.
    std::string outputImagePath = "out.png";
    if(size_t(argc) >= 3 + cgFunction.getNumArgs())
        outputImagePath = argv[3 +cgFunction.getNumArgs()-1];

    // Look at a subimages so we can process edges safely.
    std::vector<Image> inputImageViews;
    for(Image &im : inputImages)
        inputImageViews.emplace_back(
            im.getData() + padding*im.width() + padding, 
            im.width() - padding*2, im.height() - padding*2,
            im.width());

    // Perpare output image.
    Image outIm(inputImageViews[0].width(), inputImageViews[0].height(), inputImageViews[0].stride());

    // std::vector<const double*>  d; d.push_back(inputImageViews[0].getData());
    // cgFunction(d, inputImageViews[0].width(), inputImageViews[0].height(), inputImageViews[0].stride(), outView.getData());


    // Process images!
    cgFunction(inputImageViews, outIm);
    
    // Write output.
    outIm.write(outputImagePath);
    return 0;
}


