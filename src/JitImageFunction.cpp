#include "JitImageFunction.h"

namespace pixslam{

JitImageFunction::JitImageFunction(const Cell &cell, bool stdOutLogging /* = false */) : logger(stdout) {
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

void JitImageFunction::operator()(const std::vector<Image> &images, Image &out) const {
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
void JitImageFunction::operator()(const std::vector<const double *> &args, 
                int w, int h, int stride,
                double *out) const {
    generatedFunction(&args[0], w, h, stride, out); 
}

JitImageFunction::~JitImageFunction(){
    AsmJit::MemoryManager::getGlobal()->free((void*)generatedFunction); 
}

AsmJit::XmmVar JitImageFunction::functionHandler(const std::string &functionName, 
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

AsmJit::XmmVar JitImageFunction::numberHandler(const std::string &number){
    double x = std::atof(number.c_str());
    AsmJit::XmmVar xVar(compiler.newXmmVar());
    SetXmmVar(compiler, xVar, x);
    return xVar;
};

AsmJit::XmmVar JitImageFunction::symbolHandler(const std::string &name){
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

void JitImageFunction::SetXmmVar(AsmJit::X86Compiler &c, AsmJit::XmmVar &v, double d){
    // TODO: store constants in memory and load from there...
    using namespace AsmJit;
    GpVar gpreg(c.newGpVar());
    uint64_t *i = reinterpret_cast<uint64_t*>(&d);
    c.mov(gpreg, i[0]);
    c.movq(v, gpreg);
    c.unuse(gpreg);
}

void JitImageFunction::PopulateBuiltInFunctionHandlerMap(){
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



}
