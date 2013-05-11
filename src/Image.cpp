#include "Image.h"

#include "stb_image.h"
#include "stb_image_write.h"

namespace pixslam{

Image::Image(const std::string &path) : data(nullptr), w(0), h(0), ownsData(true){
    int n;
    // load image as greyscale
    unsigned char *datauc = stbi_load(path.c_str(), &w, &h, &n, 1);
    s = w;

    if(w*h > 0){
        data = new PixType[w*h];
        for(int i = 0; i < w*h; ++i)
            data[i] = datauc[i]/255.0;
    }

    stbi_image_free(datauc);
}

void Image::write(const std::string &dest) const {
    std::vector<unsigned char> datauc(w*h);

    for(int i = 0; i < w; ++i)
        for(int j = 0; j < h; ++j)
            datauc[i*w + j] = (unsigned char)(std::max(
                              std::min((*this)(i,j)*255.0, 255.0), 0.0));

    stbi_write_png(dest.c_str(), w, h, 1, &datauc[0], w);
}

}
