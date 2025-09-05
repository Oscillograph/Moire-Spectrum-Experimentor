#ifndef MOIREE_SPECTRUM_GUI_H
#define MOIREE_SPECTRUM_GUI_H

#include <include/savannah/savannah.h>

// #include <include/yaml_wrapper.h>
#include <include/logger.h>
//#include <initializer_list>

#include <external/fftw/fftw3.h> // to calculate FFT

#include <include/fwd.h>
#include <include/image.h> // for Image class

#include <thread>

namespace Savannah 
{
    struct dVec3
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };
    
    struct Signal
    {
        Signal();
        ~Signal();
        
        void Reset();
        
        int dimensions = 0; // 1 - only cols; 2 - cols and rows
        int cols = 0;
        int rows = 0;
        
        // raw data
        fftw_complex* raw_data = nullptr; // source data, is allocated only in CalculateSignal() function
        bool raw_calculated = false;
        
        // fftw data
        bool fftw_planned = false;
        fftw_complex* fftw_out_data = nullptr;
        fftw_plan fftw_plan_data;
        fftw_complex* fftw_source_row = nullptr; // updated as the need for FFTW arises
        fftw_complex* fftw_out_row = nullptr; // updated as the need for FFTW arises
        int fftw_dc_x = 0; // index of the zero-frequency corresponding to cols property
        bool fftw_x_even = false;
        int fftw_dc_y = 0; // index of the zero-frequency corresponding to rows property
        bool fftw_y_even = false;
        
        // spectrum data
        double* amplitude_spectrum = nullptr;
        double* phase_spectrum = nullptr;
        double* frequencies_axis_values = nullptr;
        
        double fftw_data_max = 0.0; // absolute value
        double fftw_data_min = 0.0; // absolute value
        bool fftw_calculated = false; // flag to escape unneeded work
        
        double amplitude_spectrum_max = 0.0; // absolute value
        double amplitude_spectrum_min = 0.0; // absolute value
    };
    
    
    struct Radar
    {
        double transmit_power = 0.0; // W
        double observation_time = 0.0; // sec
        double modulation_time = 0.0; // sec
        double deviation_frequency = 0.0; // Hz * Hz
        double starting_frequency = 0.0; // Hz
        double sampling_frequency = 0.0; // Hz
        double time_sample_width = 0.0; // sec
        dVec3 position = {0.0, 0.0, 0.0};
    };
    
    struct PointReflector
    {
        double esr = 0.0; // square meters
        double phase_shift = 0.0;
        dVec3 position = {0.0, 0.0, 0.0};
    };
    
    struct Medium
    {
        double propagation_coefficient = 1.0;
    };
    
    struct Universe
    {
        const double pi = 3.14159265358979323846;
        const double pi_half = pi/2;
        const double pi_square = pi*pi;
//        const double c = 299792458.0; // meters per second
        const double c = 3e8; // meters per second
    };
    
    struct SimulationRules
    {
        bool stop_and_go = true;
        bool no_amplitude = true;
        bool one_distance = true;
        
        bool spare_memory_on_frequencies = true;
        double beat_frequency_min = 0.0; // Hz
        double beat_frequency_max = 0.0; // Hz
    };
    
    struct Object
    {
        std::vector<PointReflector> structure = {};
        dVec3 position = {0.0, 100.0, 0.0}; // is added to all point reflectors positions
        double distance_between_elements = 0.01; // meters
        double linear_size = 100.0; // meters
        int size = 1000;
    };
    
    enum class SignalDimensionsMode
    {
        Time_Amplitude                     ,
        Time_CarrierFrequency              ,
        Time_BeatFrequency                 ,
        Time_ModulationTime                ,
        Time_ObservationTime               ,
        BeatFrequency_CarrierFrequency     ,
        BeatFrequency_CarrierFrequency_B   ,
        BeatFrequency_Height               ,
        BeatFrequency_ModulationTime       ,
        BeatFrequency_ObservationTime      ,
        BeatFrequency_Dy                   ,
    };
    
    struct GlobalUIData
    {
        bool calculation_is_going = false;
        bool calculation_is_done = false;
        bool thread_locked = false;
        std::string percentage_status = "";
        double percentage_done = 0.0;
        int todo_count = 0;
        int todo_active_index = 0;
        double signal_rows_variable_start = 0.0; // whatever the variable is, it's its value starts from here
        double signal_rows_variable_step = 0.0; // whatever the variable is, it's changed by this value each row
        double signal_rows_variable_min = 0.0; // whatever the variable is, this is its minimum value
        double signal_rows_variable_max = 0.0; // whatever the variable is, this is its maximum value
        SignalDimensionsMode signal_dimensions_mode = SignalDimensionsMode::Time_Amplitude;
    };
    GlobalUIData global_data;
    
	
	// basically, commands to process
	enum class MSAppTasks : int
	{
		Idle						,
		
        CreateImage                 ,
		ShowImage					,
		HideImage					,
		LoadImage					,
		UnloadImage					,
		
		CalculateFFT				,
		ShowFFT						,
		HideFFT						,
		SaveFFT						,
		
		ApplyColorMap				,
		ChangeColorMap				,
		
		Exit						
	};
	
	// basically, states which define what to draw on screen
	enum class MSAppMode : int
	{
		Idle						,
		
		CalculateFFT				,
		
        ShowAnalog                  ,
        SaveAnalog                  ,
		ShowFFT						,
		SaveFFT						,
		
		ChangeColorMap				, // works only for 2d spectrum images
		
		Exit
		
	};
	
	class MSApp : public App 
	{
	public:
		MSApp();
		~MSApp();
		
		// General Savannah App interface
		// void PreSetup() override;
		void SetupResources() override;
		void Logic() override;
		void GUIContent() override;
		
    // obviously, this section is not to be used by anything else in the program
	public:
		std::string m_SkillsFile = "./data/skillsDB.txt";
		std::string m_ErrorMessage = ""; // this is shown to the user
		bool m_ErrorFlag = false; // this is used to control the program flow
		bool m_FFTCalculated = false;
		
		MSAppMode m_CurrentMode = MSAppMode::Idle;
		MSAppTasks m_CurrentTask = MSAppTasks::Idle;
		std::vector<MSAppTasks> m_TaskStack = {};
		std::vector<std::string> m_ColorSchemesNamesList = {};
		int m_ColorSchemeSelected = 0;
		std::string m_ColorSchemeSelectedName = "";
		std::unordered_map<std::string, std::vector<ImVec4>> m_ColorSchemesMap = {};
		std::unordered_map<std::string, MSAppColorScheme> m_ColorSchemes = {};
		FourierSpectrumMode m_FourierSpectrumMode = FourierSpectrumMode::Amplitude;
		
		float TEXT_BASE_WIDTH = 0.0f;
		float TEXT_BASE_HEIGHT = 0.0f;
		
        Image* m_Image = nullptr;
        int m_Image_width = 1024;
        int m_Image_height = 1024;
		
        // data
        Universe universe;
        Medium medium;
        Radar radar;
        Object object;
        Signal signal;
        SimulationRules rules;
        
		// color features
		float lowerFFTWLevel = 0.0f;
		float upperFFTWlevel = 0.5f;
		
		void NewTask(MSAppTasks task);
		
		void ShowMainMenu();
		void ShowContent();
	};
    
	void CalculateSignal(Radar* radar, Object* object, Medium* medium, Universe* universe, SimulationRules* rules, Signal* signal);
    void CalculateSpectrum(Signal* signal, Radar* radar);
    void CalculateThread(MSApp* app, Radar* radar, Object* object, Medium* medium, Universe* universe, SimulationRules* rules, Signal* signal);
    
//    void CalculateDistances(Radar* radar, Object* object, SimulationRules* rules, double* distances = nullptr);
    
    Image* CreateImage(MSApp* app, Image* image, Signal& signal, ImVec2 image_size_desired = {640, 480});
}

#endif
