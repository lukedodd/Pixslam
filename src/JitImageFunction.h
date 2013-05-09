#pragma once

#include <functional>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <map>

#include "asmjit/asmjit.h"


#include "Image.h"
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

class JitImageFunction : public Visitor<AsmJit::XmmVar>{
public:
    JitImageFunction(const Cell &cell, bool stdOutLogging = false);

    // Call JIT function.
    void operator()(const std::vector<Image> &images, Image &out) const;

    // "Lower level" call for image data from other sources (e.g. opencv)
    void operator()(const std::vector<const double *> &args, 
                    int w, int h, int stride,
                    double *out) const;

    virtual size_t getNumArgs() const {return argNameToIndex.size();}

    virtual ~JitImageFunction();

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
            const std::vector<AsmJit::XmmVar> &args);

    virtual AsmJit::XmmVar numberHandler(const std::string &number);

    virtual AsmJit::XmmVar symbolHandler(const std::string &name);

    void SetXmmVar(AsmJit::X86Compiler &c, AsmJit::XmmVar &v, double d);

    void PopulateBuiltInFunctionHandlerMap();

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
