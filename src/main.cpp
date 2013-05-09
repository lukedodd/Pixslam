#include <iostream>
#include <sstream>
#include <fstream>


#include "Image.h"
#include "Parser.h"
#include "Compiler.h"

using namespace pixslam;


void logCommandLine(int argc, char *argv[], const std::string &filePrefix){
    #ifdef _WIN32
        std::string fileName = filePrefix + ".bat";
    #else
        std::string fileName = filePrefix + ".sh";
    #endif
    std::ofstream out(fileName, std::ios::out);

    for(int i = 0; i < argc; ++i)
        out << argv[i] << " ";

    out << std::endl;
}

int main (int argc, char *argsRaw[])
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

    std::vector<std::string> argv;
    std::vector<std::string> options;

    // separate options from file/code arguments
    for(int i = 0; i < argc; ++i){
        char *s = argsRaw[i];
        if(*s != '-')
            argv.emplace_back(s);
        else
            options.emplace_back(s);
    }

    // parse command line arguments
    bool logAsm = false;
    bool logCommand = false;
    for(auto s : options){
        if(s == "--logAsm")
            logAsm = true;
        else if(s == "--logCommand")
            logCommand = true;
        else{
            std::cerr << "Unrecognised command line switch: " << s << std::endl;
            return 1;
        }
            
    }

    // See if first arg is a file and read code from it.
    // We infer file by checking that first char is not a '(' (cheeky but it works!)
    // Otherwise we interperate the argument as code directly.
    std::string codeString;
    if(argv[1].size() > 0 && argv[1][0] != '('){
        std::ifstream ifs(argv[1]);
        if(ifs){
            std::stringstream buffer;
            buffer << ifs.rdbuf();
            codeString = buffer.str();
        }else{
            std::cout << "Could not find file " << argv[1] << std::endl;
            return 1;
        }
    }else{ 
        codeString = argv[1];
    }

    // Generate code.
    Cell code = cellFromString(codeString);
    Compiler cgFunction(code, logAsm);

    // Read in input images specified by arguments.
    int padding = 5;
    std::vector<Image> inputImages;
    for(size_t i = 0; i < cgFunction.getNumArgs(); ++i){
        Image im(argv[2+i]);

        if(im.width()*im.height() == 0){
            std::cout << "Failed to load image " << argv[2+i] << std::endl;
            return 1;
        }

        inputImages.emplace_back(im, padding, padding);
    }

    // Remaining arg, if preset is our output destination.
    std::string outputImagePath = "out.png";
    if(size_t(argc) >= 3 + cgFunction.getNumArgs())
        outputImagePath = argv[3 +cgFunction.getNumArgs()-1];

    // log command line if requested (useful for easy to understand examples directory)
    if(logCommand) logCommandLine(argc, argsRaw, outputImagePath);

    // Look at a subimages so we can process edges safely.
    std::vector<Image> inputImageViews;
    for(Image &im : inputImages)
        inputImageViews.emplace_back(
            im.getData() + padding*im.width() + padding, 
            im.width() - padding*2, im.height() - padding*2,
            im.width());

    // Perpare output image.
    Image outIm(inputImageViews[0].width(), inputImageViews[0].height(), inputImageViews[0].stride());

    // Process images!
    cgFunction(inputImageViews, outIm);
    
    // Write output.
    outIm.write(outputImagePath);
    return 0;
}


