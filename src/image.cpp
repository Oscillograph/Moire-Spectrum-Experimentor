#include <include/image.h>

//#define STB_IMAGE_IMPLEMENTATION
#include <external/stb_image.h>
#include <external/stb_image_write.h>
#include <filesystem>

#include <include/savannah/platforms/opengl/opengl_texture.h> // for OpenGLTexture2D

namespace Savannah
{
    Image::Image()
    {}
    
    Image::~Image()
    {
        Unload();
        
        delete texture;
        texture = nullptr;
        
        delete originalTexture;
        originalTexture = nullptr;
        
        delete[] fftw_data;
        fftw_data = nullptr;
        
        delete[] pixels;
        pixels = nullptr;
    }
    
    void Image::Create(int w, int h, void* data)
    {
        width = w;
        height = h;
        channelsOriginal = 4;
        channelsDesired = 4;
        
        if (data != nullptr)
        {
            // prepare image from data
            pixels = (unsigned char*)data;
            
            UpdateTexture();
//            UpdateOriginalTexture();
            CONSOLE_LOG("Texture created");
        }
    }
    
    void Image::Load(const std::string& filepath)
    {
        // load pixels
        int w, h, c;
        pixels = stbi_load(filepath.c_str(), &w, &h, &c, channelsDesired);
        if (!pixels)
        {
            CONSOLE_LOG("stbi_load returned NULL. File format is probably not suppored.");
            return;
        } else {
            width = w;
            height = h;
            channelsOriginal = c;
        }
        
        aspectRatio = (float)height / width;
        CONSOLE_LOG("Loaded Image constraints (w, h, c, a/r): ", width, ", ", height, ", ", channelsOriginal, ", ", aspectRatio);
        CONSOLE_LOG("Desired channels: ", channelsDesired);
        
        // reinitialize values
        fftw_data_max = 0.0;
        fftw_data_min = 1.0;
        magnitudeOrderPolishCoefficient = 1.0;
        
        // load fftw data from pixels
        CONSOLE_LOG("Generating FFTW data");
        delete[] fftw_data;
        fftw_data = new fftw_complex[width * height];
        unsigned char* pixel = pixels;
        uint32_t p, r, g, b, a;
        p = 256 * 256 * 256 * 256 - 1;
        r = 256 * 256 * 256;
        g = 256 * 256;
        b = 256;
        a = 1;
        for (int x = 0; x < width; ++x)
        {
            int xIndex = x*height; 
            for (int yIndex = 0; yIndex < height; ++yIndex)
            {
                if (channelsOriginal == 4)
                {
                    uint32_t value = *(pixels + (xIndex + yIndex) * channelsOriginal) * r + 
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 1) * g +
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 2) * b +
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 3) * a;
                    
                    fftw_data[xIndex + yIndex][0] = (float)value / p;
                    
                    if (colors.find(value) != colors.end())
                    {
                        colors[value] += 1;
                    } else {
                        colors[value] = 1;
                    }
                }
                if (channelsOriginal == 3)
                {
                    uint32_t value = *(pixels + (xIndex + yIndex) * channelsOriginal) * g + 
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 1) * b +
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 2) * a;
                    
                    fftw_data[xIndex + yIndex][0] = (float)value / r;
                    
                    if (colors.find(value) != colors.end())
                    {
                        colors[value] += 1;
                    } else {
                        colors[value] = 1;
                    }
                }
            }
        }
        CONSOLE_LOG("FFTW data generation done");
        
        // load texture to display
        UpdateTexture();
        UpdateOriginalTexture();
        ready = true;
        CONSOLE_LOG("Texture updated");
    }
    
    void Image::UpdatePixels(MSAppColorScheme* colorScheme, FourierSpectrumMode mode)
    {
        if (!pixels)
        {
            return;
        }
        if (!colorScheme)
        {
            CONSOLE_RED("Image::UpdatePixels: color scheme is nullptr");
            return;
        } else {
            CONSOLE_LOG("Chosen scheme: ", colorScheme->name);
        }
        magnitudeOrderPolishCoefficient = brightnessCoefficient*sqrt((fftw_data_max / fftw_data_min));
        magnitudeOrderPolishCoefficient = (magnitudeOrderPolishCoefficient == INFINITY) ? 100.0 : magnitudeOrderPolishCoefficient;		
        CONSOLE_LOG("FFTW Data. Max: ", fftw_data_max, "; Min: ", fftw_data_min);
        CONSOLE_LOG("Magnitude Order: ", magnitudeOrderPolishCoefficient);
        
        unsigned char* pixel = pixels;
        for (int x = 0; x < width; ++x)
        {
            int xIndex = x*height;
            int red, green, blue = 0;
            
            for (int yIndex = 0; yIndex < height; ++yIndex)
            {
                double value = 0.0;
                
                switch (mode)
                {
                case FourierSpectrumMode::Amplitude:
                {
                    value = sqrt(fftw_data[xIndex + yIndex][0] * fftw_data[xIndex + yIndex][0] + fftw_data[xIndex + yIndex][1] * fftw_data[xIndex + yIndex][1]);
                    value = (value * magnitudeOrderPolishCoefficient < 1.0) ? value * magnitudeOrderPolishCoefficient : 1.0;
                }
                    break;
                case FourierSpectrumMode::Phase:
                    {
                        if (fftw_data[xIndex + yIndex][0] == 0)
                        {
                            if (fftw_data[xIndex + yIndex][1] >= 0)
                            {
                                value = pi_half;
                            } else {
                                value = -pi_half;
                            }
                        } else {
                            value = atan(fftw_data[xIndex + yIndex][1] / fftw_data[xIndex + yIndex][0]);
                        }
                        
                        value = (value + pi_half) / pi;
                    }
                    break;
                case FourierSpectrumMode::Real:
                    {
                        value = fftw_data[xIndex + yIndex][0];
                        value = (value * magnitudeOrderPolishCoefficient < 1.0) ? value * magnitudeOrderPolishCoefficient : 1.0;
                    }
                    break;
                case FourierSpectrumMode::Imaginary:
                    {
                        value = fftw_data[xIndex + yIndex][1];
                        value = (value * magnitudeOrderPolishCoefficient < 1.0) ? value * magnitudeOrderPolishCoefficient : 1.0;
                    }
                    break;
                }
                
                colorScheme->GetColor(value, fftw_data_min, fftw_data_max, &red, &green, &blue);
                
                *(pixels + (xIndex + yIndex) * channelsOriginal) = red;
                *(pixels + (xIndex + yIndex) * channelsOriginal + 1) = green;
                *(pixels + (xIndex + yIndex) * channelsOriginal + 2) = blue;
                if (channelsOriginal == 4)
                    *(pixels + (xIndex + yIndex) * channelsOriginal + 3) = 255; // alpha
            }
        }
    }
    
    void Image::UpdateTexture()
    {
        if (!pixels)
        {
            return;
        }
        delete texture;
        CONSOLE_LOG("Creating texture for Image object");
        texture = new OpenGLTexture2D(pixels, width, height, channelsOriginal);
    }
    
    void Image::UpdateOriginalTexture()
    {
        if (!pixels)
        {
            return;
        }
        delete originalTexture;
        CONSOLE_LOG("Creating \"original\" texture for Image object");
        originalTexture = new OpenGLTexture2D(pixels, width, height, channelsOriginal);
    }
    
    void Image::NormalizeFFT()
    {
        // update fftw data absolute min and max metrics
        fftw_data_max = 0.0; // absolute value
        fftw_data_min = 1.0; // absolute value
        for (int x = 0; x < width; ++x)
        {
            for (int y = 0; y < height; ++y)
            {
                double value;
                
                value = fabs(fftw_data[x*height + y][0]);
                if ((value > fftw_data_max) && value != 1.0)
                    fftw_data_max = value;
                if ((value < fftw_data_min) && value != 0.0)
                    fftw_data_min = value;
                
                value = fabs(fftw_data[x*height + y][1]);
                if ((value > fftw_data_max) && value != 1.0)
                    fftw_data_max = value;
                if ((value < fftw_data_min) && value != 0.0)
                    fftw_data_min = value;
            }
        }
        
        // normalize
        for (int x = 0; x < width; ++x)
        {
            for (int y = 0; y < height; ++y)
            {
                fftw_data[x*height + y][0] = fftw_data[x*height + y][0] / fftw_data_max;
                fftw_data[x*height + y][1] = fftw_data[x*height + y][1] / fftw_data_max;
            }
        }
        
        fftw_data_max = 1.0;
    }
    
    void Image::Unload()
    {
        if (pixels)
        {
            stbi_image_free(pixels);
            pixels = nullptr;
            colors.clear();
            colorsSorted.clear();
//			width = 0;
//			height = 0;
//			channelsOriginal = 0;
//			aspectRatio = 1.0;
        }
    }
    
    void Image::Save(const std::string& path)
    {
        CONSOLE_LOG("stbi call to write png");
        int stride_in_bytes = width * channelsOriginal;
        stbi_write_png(path.c_str(), width, height, channelsOriginal, pixels, stride_in_bytes);
    }
    
    void MSAppColorScheme::GetColor(double value, double min, double max, int* r, int* g, int* b)
    {
        int left, right;
//		CONSOLE_LOG(colors.size());
        
        if (value <= min)
        {
            left = 0;
            right = 0;
        } else if (value >= max) {
            left = colors.size() - 1;
            right = colors.size() - 1;
        } else {
            // "move" value across the color map into appropriate region
            value = value * (colors.size() - 1);
            left = (int)floor(value);
            right = left + 1;
            value = value - (double)left;
        }
        
        *r = (colors[right].x - colors[left].x) * value + colors[left].x;
        *g = (colors[right].y - colors[left].y) * value + colors[left].y;
        *b = (colors[right].z - colors[left].z) * value + colors[left].z;
        
        return;
    }
}
