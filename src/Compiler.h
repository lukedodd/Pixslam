#pragma once

#include <functional>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <map>

#include "asmjit/asmjit.h"

#include "Parser.h"

namespace pixslam{

template <typename EvalReturn> class Visitor{
public:
    typedef const double * const * Arguments;
    typedef void (*FuncPtrType)(Arguments, size_t, size_t, size_t, double *);


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

class Compiler : public Visitor<AsmJit::XmmVar>{
public:
    Compiler(const Cell &cell, bool stdOutLogging = false) : logger(stdout) {

        if(stdOutLogging)
            compiler.setLogger(&logger);
       
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

    virtual ~Compiler(){
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

        // std::cout << argNameToIndex.size() << std::endl;
        // std::cout << argNameToIndex["A"] << std::endl;
        // std::cout << "functionName: " << functionName << std::endl;

        // try builtin function lookup first
        auto it = functionHandlerMap.find(functionName);
        if(it != functionHandlerMap.end()){
            return it->second(args);
        }

        // std::cout << "Not a function." << std::endl;

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

            compiler.cmpsd(ret, args[1], 1); // LT
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap[">"] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 1); // LT
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap["<="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[0]);

            compiler.cmpsd(ret, args[1], 2); // LEQ
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap[">="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 2); // LEQ
            compiler.andpd(ret, one);
            return ret;
        };
            
        functionHandlerMap["=="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 0); // EQ
            compiler.andpd(ret, one);
            return ret;
        };

        functionHandlerMap["!="] = [&](const std::vector<XmmVar> &args) -> XmmVar{
            XmmVar ret = compiler.newXmmVar();
            compiler.movq(ret, args[1]);

            compiler.cmpsd(ret, args[0], 4); // NEQ
            compiler.andpd(ret, one);
            return ret;
        };
    }

private:

    typedef std::function<AsmJit::XmmVar (const std::vector<AsmJit::XmmVar> &)> 
        BuiltInFunctionHandler;
    std::map<std::string, BuiltInFunctionHandler> functionHandlerMap;

    AsmJit::X86Compiler compiler;
    std::map<std::string, size_t> argNameToIndex;

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

};

}
