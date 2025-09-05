#ifndef MSAPP_IMAGE_H
#define MSAPP_IMAGE_H

#include <include/fwd.h>
#include <external/imgui/imgui.h> // for ImVec4
#include <include/logger.h>
#include <external/fftw/fftw3.h> // to calculate FFT

#include <unordered_map>

namespace Savannah
{
    class OpenGLTexture2D;
    
    struct MSAppColorScheme
    {
        std::string name = "";
        std::vector<ImVec4> colors = {};
        
        void GetColor(double value, double min, double max, int* r, int* g, int* b);
    };
    
    class Image
    {
    public:
        Image();
        ~Image();
        
        void Create(int w, int h, void* data = nullptr);
        void Load(const std::string& filepath);
        void UpdatePixels(MSAppColorScheme* colorScheme, FourierSpectrumMode mode = FourierSpectrumMode::Amplitude); // align pixels data with fftw data
        void UpdateTexture(); // if pixels were manipulated
        void UpdateOriginalTexture(); // if an image was loaded
        void NormalizeFFT();
        void Unload();
        void Save(const std::string& path);
        
        int width = 0;
        int height = 0;
        float aspectRatio = 1.0;
        int channelsOriginal = 0; // bit per pixel
        int channelsDesired = 0;
        double brightnessCoefficient = 2.0; // can be set up by user
        double magnitudeOrderPolishCoefficient = 1.0; // calculated dynamically by the app
        std::string path = "";
        std::unordered_map<uint32_t, uint32_t> colors = {};
        std::vector<std::pair<uint32_t, uint32_t>> colorsSorted = {};
        
        fftw_complex* fftw_data = nullptr;
        double fftw_data_max = 0.0; // absolute value
        double fftw_data_min = 0.0; // absolute value
        unsigned char* pixels = nullptr;
        OpenGLTexture2D* texture = nullptr;
        OpenGLTexture2D* originalTexture = nullptr;
        bool ready = false;
    };
}

#endif
