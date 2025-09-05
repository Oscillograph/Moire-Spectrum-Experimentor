#include <include/MoireSpectrum_gui.h>

//#include <src/savannah/platforms/opengl/opengl_texture.cpp>

namespace Savannah 
{	
    Signal::Signal()
    {
    }
    
    Signal::~Signal()
    {
        Reset();
    }
    
    void Signal::Reset()
    {
        delete[] raw_data;
        raw_data = nullptr;
        this->raw_calculated = false;
        
        fftw_destroy_plan(fftw_plan_data);
        fftw_free(fftw_out_data);
        fftw_out_data = nullptr;
        this->fftw_planned = false;
        
        delete[] amplitude_spectrum;
        amplitude_spectrum = nullptr;
        delete[] phase_spectrum;
        phase_spectrum = nullptr;
        delete[] frequencies_axis_values;
        frequencies_axis_values = nullptr;
        this->fftw_calculated = false;
    }
    
    void CalculateSignal(Radar* radar, 
                         Object* object, 
                         Medium* medium, 
                         Universe* universe, 
                         SimulationRules* rules, 
                         Signal* signal)
    {
        CONSOLE_LOG("CalculateSignal - start");
        signal->raw_calculated = false;
        if (signal->dimensions > 0) // making sure the signal is initialized at all
        {
            // spare resouces
            if (rules->spare_memory_on_frequencies)
            {
                rules->beat_frequency_min = 0.0;
                
                // this memory saving technique is based on the assumption that the object is basically a line perpendicular to the X axis
                double distance_max = 0.0;
                double distance = 0.0;
                for (int i = 0; i < object->size; ++i)
                {
                    distance = sqrt(pow(radar->position.x - (object->position.x + object->structure[object->size - 1].position.x), 2) +
                                    pow(radar->position.y - (object->position.y + object->structure[object->size - 1].position.y), 2) +
                                    pow(radar->position.z - (object->position.z + object->structure[object->size - 1].position.z), 2));
                    if (distance > distance_max)
                    {
                        distance_max = distance;
                    }
                }
                double tk_max = 2.1*distance_max / universe->c;
                rules->beat_frequency_max = radar->deviation_frequency * tk_max;
                radar->sampling_frequency = 2 * rules->beat_frequency_max;
                radar->time_sample_width = 1/radar->sampling_frequency;
                CONSOLE_LOG("Sparing memory mode ON");
                CONSOLE_LOG("tk_max: ", tk_max, " sec");
                CONSOLE_LOG("beat frequency max: ", rules->beat_frequency_max, " Hz");
                CONSOLE_LOG("sampling frequency: ", radar->sampling_frequency, " Hz");
            } else {
                // still spare memory, but not that hard
                rules->beat_frequency_max = radar->deviation_frequency * radar->observation_time;
                radar->sampling_frequency = 2 * rules->beat_frequency_max;
                radar->time_sample_width = 1/radar->sampling_frequency;
                CONSOLE_LOG("Sparing memory mode OFF");
                CONSOLE_LOG("beat frequency max: ", rules->beat_frequency_max, " Hz");
                CONSOLE_LOG("sampling frequency: ", radar->sampling_frequency, " Hz");
            }
            
            // Important! cols might become so big that signal data exceed memory available
            int cols = (int)(round(radar->observation_time / radar->time_sample_width));
            if (signal->cols != cols)
            {
                if (signal->fftw_planned)
                {
                    signal->fftw_planned = false;
                    fftw_destroy_plan(signal->fftw_plan_data);
                    fftw_free(signal->fftw_out_data);
                    signal->fftw_out_data = nullptr;
                }
                signal->cols = cols;
            }
            
            if (global_data.signal_dimensions_mode == SignalDimensionsMode::Time_Amplitude)
            {
                signal->rows = 1;
            } else {
                if ((global_data.signal_rows_variable_max != 0.0) && 
                    (global_data.signal_rows_variable_step != 0.0))
                {
                    signal->rows = (int)floor(round((global_data.signal_rows_variable_max - global_data.signal_rows_variable_min) / global_data.signal_rows_variable_step));
                    CONSOLE_LOG("signal rows to be calculated: ", signal->rows);
                } else {
                    signal->rows = 1;
                }
            }
            
            CONSOLE_LOG("total samples: ", signal->cols, " x ", signal->rows);
            CONSOLE_LOG("total memory usage: ", (int)(round(signal->cols * sizeof(fftw_complex) * signal->rows / 1024)), " KB");
            
            if (signal->raw_data != nullptr)
            {
                delete[] signal->raw_data;
            }
            int fftw_elements_total = signal->rows * signal->cols;
            signal->raw_data = new fftw_complex[fftw_elements_total];
            for (size_t i = 0; i < fftw_elements_total; ++i)
            {
                signal->raw_data[i][0] = 0.0;
                signal->raw_data[i][1] = 0.0;
            }
            signal->fftw_source_row = signal->raw_data; // temporary measure to properly init fftw_plans
            // the idea is to navigate on "raw_data" using "fftw_source_row += cols"
            // it based on the assumption that FFTW library stores a pointer to the source data, not its value
            // same thing goes with fftw_out_row
//            CONSOLE_LOG("zero: ", global_data.zero, "; zero2: [", global_data.zero2[0], ", ", global_data.zero2[1], "]");
            if (!signal->fftw_planned)
            {
                signal->fftw_out_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * signal->rows * signal->cols);
                
                for (size_t i = 0; i < fftw_elements_total; ++i)
                {
                    signal->fftw_out_data[i][0] = 0.0;
                    signal->fftw_out_data[i][1] = 0.0;
                }
                signal->fftw_out_row = signal->fftw_out_data;
                
                signal->fftw_plan_data = fftw_plan_dft_1d(cols, signal->fftw_source_row, signal->fftw_out_row, FFTW_FORWARD, FFTW_ESTIMATE);
                signal->fftw_planned = true;
            }
            
            double* time_samples = new double[signal->cols];
            double* distance = new double[signal->cols];
            double* amplitude = new double[signal->cols];

            for (size_t m = 0; m < signal->cols * signal->rows; ++m)
            {
                if (m < signal->cols)
                {
                    amplitude[m] = 0.0;
                    distance[m] = 0.0;
                    time_samples[m] = m * radar->time_sample_width;
                }
                
                signal->raw_data[m][0] = 0.0;
                signal->raw_data[m][1] = 0.0;
                signal->fftw_out_data[m][0] = 0.0;
                signal->fftw_out_data[m][1] = 0.0;
            }
            
            CONSOLE_LOG("object consists of ", object->structure.size(), " elements");
            for (int object_index = 0; object_index < object->structure.size(); ++object_index)
            {
                // less calculations
                if (global_data.signal_dimensions_mode != SignalDimensionsMode::BeatFrequency_Height)
                if (rules->one_distance)
                {
                    double one_distance = 2*sqrt(
                                               pow((radar->position.x - (object->position.x + object->structure[object_index].position.x)), 2) +
                                               pow((radar->position.y - (object->position.y + object->structure[object_index].position.y)), 2) + 
                                               pow((radar->position.z - (object->position.z + object->structure[object_index].position.z)), 2)
                                               );
                    for (int j = 0; j < signal->cols; ++j)
                    { 
                        distance[j] = one_distance;
                    }
                } else {
                    // different distances between transmit antenna and object, and receiver antenna and object
                    double distance1 = 0.0;
                    double distance2 = 0.0;
                    for (int j = 0; j < signal->cols; ++j)
                    {
                        distance[j] = distance1 + distance2;
                    }
                }
                
                // makes sense for rows > 1
                if (rules->stop_and_go)
                {}
            
                // less calculations
                if (rules->no_amplitude)
                {
                    for (int j = 0; j < signal->cols; ++j)
                    {
                        amplitude[j] = 1.0;
                    }
                } else {
                    for (int j = 0; j < signal->cols; ++j)
                    {
                        amplitude[j] = 0.0; // TODO: make use of distances
                    }
                }
                
                switch (global_data.signal_dimensions_mode)
                {
                    case SignalDimensionsMode::Time_Amplitude:
                        {
//                            CONSOLE_LOG("Time_Amplitude domain calculation for obj=", object_index);
                            double beat_frequency = 0.0;
                    
                            for (int row_index = 0; row_index < signal->rows; ++row_index)
                            {
                                signal->fftw_source_row = signal->raw_data + row_index * signal->cols;
                                for (int j = 0; j < signal->cols; ++j)
                                {
                                    double tk = distance[j] / universe->c;
                                    beat_frequency = 0.5 * radar->deviation_frequency * tk;
        //                        CONSOLE_LOG("tk = ", tk);
        //                        CONSOLE_LOG("beta*tk*t = ", radar->deviation_frequency * tk * time_samples[j]);
                                    
                                    // real part
                                    double a1 = - 2.0 * beat_frequency * time_samples[j];
                                    double a2 = - radar->starting_frequency * tk;
                                    double a3 = beat_frequency * tk;
                                    double argument = 2 * pi * (a1 + a2 + a3);
                                    double amp = 0.5 * amplitude[j] * medium->propagation_coefficient;
                                    double result = amp * cos(argument);
                                    signal->fftw_source_row[j][0] = signal->fftw_source_row[j][0] + result;
                                    
                                    // imaginary part
                                    signal->fftw_source_row[j][1] += 0.0;
                                }
                            }
                        }
                        break;
                    case SignalDimensionsMode::BeatFrequency_CarrierFrequency:
                        {
//                            CONSOLE_LOG("BeatFrequency_CarrierFrequency domain calculation");
                            double carrier_frequency = global_data.signal_rows_variable_min;
                            double beat_frequency = 0.0;
                            
                            for (int row_index = 0; row_index < signal->rows; ++row_index)
                            {
                                int fftw_source_row_index_start = row_index * signal->cols;
                                carrier_frequency = global_data.signal_rows_variable_min + row_index * global_data.signal_rows_variable_step;
                                signal->fftw_source_row = signal->raw_data + fftw_source_row_index_start;
                                
                                for (int j = 0; j < signal->cols; ++j)
                                {
                                    double tk = distance[j] / universe->c;
                                    beat_frequency = 0.5 * radar->deviation_frequency * tk;

                                    // real part
                                    double a1 = - 2.0 * beat_frequency * time_samples[j];
                                    double a2 = - carrier_frequency * tk;
                                    double a3 = beat_frequency * tk;
                                    double argument = 2 * pi * (a1 + a2 + a3);
                                    double amp = 0.5 * amplitude[j] * medium->propagation_coefficient;
                                    double result = amp * cos(argument);
                                    signal->fftw_source_row[j][0] += result;

                                    // imaginary part
//                                    signal->fftw_source_row[j][1] += 0.0;
                                }
                            }
                        }
                        break;
                    case SignalDimensionsMode::BeatFrequency_Height:
                        {
//                            CONSOLE_LOG("BeatFrequency_CarrierFrequency domain calculation");
                            double height = global_data.signal_rows_variable_min;
                            double beat_frequency = 0.0;
                            
                            for (int row_index = 0; row_index < signal->rows; ++row_index)
                            {
                                int fftw_source_row_index_start = row_index * signal->cols;
                                signal->fftw_source_row = signal->raw_data + fftw_source_row_index_start;
                                height = global_data.signal_rows_variable_min + row_index * global_data.signal_rows_variable_step;
                                
                                // calculate distances to objects
                                if (rules->one_distance)
                                {
                                    double one_distance = 2*sqrt(
                                                                 pow((radar->position.x - (object->position.x + object->structure[object_index].position.x)), 2) +
                                                                 pow((radar->position.y - (object->position.y + object->structure[object_index].position.y)), 2) + 
                                                                 pow((height - (object->position.z + object->structure[object_index].position.z)), 2)
                                                                 );
                                    for (int j = 0; j < signal->cols; ++j)
                                    { 
                                        distance[j] = one_distance;
                                    }
                                } else {
                                    // different distances between transmit antenna and object, and receiver antenna and object
                                    double distance1 = 0.0;
                                    double distance2 = 0.0;
                                    for (int j = 0; j < signal->cols; ++j)
                                    {
                                        distance[j] = distance1 + distance2;
                                    }
                                }
                                
                                // calculate signal
                                for (int j = 0; j < signal->cols; ++j)
                                {
                                    double tk = distance[j] / universe->c;
                                    beat_frequency = 0.5 * radar->deviation_frequency * tk;

                                    // real part
                                    double a1 = - 2.0 * beat_frequency * time_samples[j];
                                    double a2 = - radar->starting_frequency * tk;
                                    double a3 = beat_frequency * tk;
                                    double argument = 2 * pi * (a1 + a2 + a3);
                                    double amp = 0.5 * amplitude[j] * medium->propagation_coefficient;
                                    double result = amp * cos(argument);
                                    signal->fftw_source_row[j][0] += result;

                                    // imaginary part
//                                    signal->fftw_source_row[j][1] += 0.0;
                                }
                            }                            
                        }
                        break;
                }
                
                // tell how much work is done
                global_data.percentage_done = 100 * (double)(object_index + 1) / object->structure.size();
            }
            delete[] distance;
            delete[] amplitude;
            delete[] time_samples;
            
            // update min & max values
            CONSOLE_LOG("Update signal min & max values");
            signal->fftw_data_min = 0.0;
            signal->fftw_data_max = 0.0;
            for (size_t i = 0; i < signal->rows; ++i)
            {
                signal->fftw_source_row = signal->raw_data + (i * signal->cols);
                for (size_t j = 0; j < signal->cols; ++j)
                {
                    if (std::abs(signal->fftw_source_row[j][0]) > std::abs(signal->fftw_data_max))
                    {
                        signal->fftw_data_max = signal->fftw_source_row[j][0];
                    }
                    
                    if (std::abs(signal->fftw_source_row[j][0]) < std::abs(signal->fftw_data_min))
                    {
                        signal->fftw_data_min = signal->fftw_source_row[j][0];
                    }
                }
            }
            
            CONSOLE_LOG("Signal calculation complete");
            signal->raw_calculated = true;
            return;
        }
        
        signal->raw_calculated = false;
    }
    
    void CalculateSpectrum(Signal* signal, Radar* radar)
    {        
        CONSOLE_LOG("CalculateSpectrum - start");
        signal->fftw_calculated = false;
        if ((nullptr != signal->amplitude_spectrum) && 
            (nullptr != signal->phase_spectrum) && 
            (nullptr != signal->frequencies_axis_values))
        {
            // clear traces of pevious calculations
            CONSOLE_LOG("CalculateSpectrum: clear traces of pevious calculations");
            delete[] signal->amplitude_spectrum;
            signal->amplitude_spectrum = nullptr;
            delete[] signal->phase_spectrum;
            signal->phase_spectrum = nullptr;
            delete[] signal->frequencies_axis_values;
            signal->frequencies_axis_values = nullptr;
            CONSOLE_LOG("Done");
        }
        
        // allocate memory
        CONSOLE_LOG("Allocating memory");
        signal->amplitude_spectrum = new double[signal->rows * signal->cols];
        signal->phase_spectrum = new double[signal->rows * signal->cols];
        signal->frequencies_axis_values = new double[signal->cols];
  
        CONSOLE_LOG("Initializing memory");
        int total_elements = signal->rows * signal->cols;
        for (int j = 0; j < total_elements; ++j)
        {
            signal->amplitude_spectrum[j] = 0.0;
            signal->phase_spectrum[j] = 0.0;
            if (j < signal->cols)
            {
                signal->frequencies_axis_values[j] = 0.0;
            }
        }
//        double zero1d = 0.0;
//        memcpy(signal->amplitude_spectrum, &zero1d, signal->rows * signal->cols * sizeof(double));
//        memcpy(signal->phase_spectrum, &zero1d, signal->rows * signal->cols * sizeof(double));
//        memcpy(signal->frequencies_axis_values, &zero1d, signal->cols * sizeof(double));
 
        // calculate spectrums
        for (int row_index = 0; row_index < signal->rows; ++row_index)
        {
            signal->fftw_source_row = signal->raw_data + (row_index * signal->cols);
            signal->fftw_out_row = signal->fftw_out_data + (row_index * signal->cols);
//            CONSOLE_CYAN("row ", row_index, ": fftw_source_row = ", signal->fftw_source_row);
//            CONSOLE_GREEN("row ", row_index, ": fftw_out_row = ", signal->fftw_out_row);
            
//            fftw_destroy_plan(signal->fftw_plan_data);
//            signal->fftw_plan_data = fftw_plan_dft_1d(signal->cols, 
//                                                      signal->fftw_source_row, 
//                                                      signal->fftw_out_row, 
//                                                      FFTW_FORWARD, 
//                                                      FFTW_ESTIMATE);
            
            fftw_execute_dft(signal->fftw_plan_data, 
                             signal->fftw_source_row, 
                             signal->fftw_out_row);
            
            // tell how much work is done
            global_data.percentage_done = 50 * (double)(row_index + 1) / signal->rows;
        }
//        CONSOLE_LOG("CalculateSpectrum - almost-end");
//        CONSOLE_LOG("Raw data cols", signal->cols);
        double frequency_axis_step = radar->sampling_frequency / signal->cols;
        signal->fftw_x_even = ((signal->cols % 2) == 0);
        signal->fftw_dc_x = (signal->fftw_x_even) ? (signal->cols / 2) : (signal->cols / 2 + 1);
        
        for (int row_index = 0; row_index < signal->rows; ++row_index)
        {
            int row_offset = row_index * signal->cols;
            signal->fftw_out_row = signal->fftw_out_data + row_offset;
            for (int i = 0; i < signal->cols; ++i)
            {
                signal->amplitude_spectrum[row_offset + i] = sqrt(pow(signal->fftw_out_row[i][0], 2) + pow(signal->fftw_out_row[i][1], 2));

                if (signal->fftw_out_row[i][0] == 0.0)
                {
                    signal->phase_spectrum[row_offset + i] = (signal->fftw_out_row[i][1] > 0.0) ? 1.0 : -1.0;
                } else {
                    signal->phase_spectrum[row_offset + i] = signal->fftw_out_row[i][1] / signal->fftw_out_row[i][0];
                }

                signal->frequencies_axis_values[i] = i * frequency_axis_step;
            }
            
            // fftshift amplitude spectrum
//            double frequencies_axis_max = signal->frequencies_axis_values[signal->cols - 1];
//            double* buffer = new double[signal->fftw_dc_x];
//            for (int i = 0; i < signal->fftw_dc_x; ++i)
//            {
//                buffer[i] = signal->amplitude_spectrum[row_offset + i];
//            }
//            if (!signal->fftw_x_even)
//            {
//                for (int i = 0; i < signal->fftw_dc_x; ++i)
//                {
//                    signal->amplitude_spectrum[row_offset + i] = 
//                                    signal->amplitude_spectrum[row_offset + signal->fftw_dc_x + i];
//                }
//                for (int i = 0; i < signal->fftw_dc_x; ++i)
//                {
//                    signal->amplitude_spectrum[row_offset + signal->fftw_dc_x + i] = buffer[i]; 
//                }
//            } else {
//                for (int i = 0; i < signal->fftw_dc_x - 1; ++i)
//                {
//                    signal->amplitude_spectrum[row_offset + i] = 
//                                    signal->amplitude_spectrum[row_offset + signal->fftw_dc_x + i];
//                }
//                for (int i = 0; i < signal->fftw_dc_x; ++i)
//                {
//                    signal->amplitude_spectrum[row_offset + signal->fftw_dc_x - 1 + i] = buffer[i]; 
//                }
//            }
//            delete[] buffer;
            
            // tell how much work is done
            global_data.percentage_done = 50.0 + (50 * (double)(row_index + 1) / signal->rows);
        }
        
        CONSOLE_LOG("Spectrum calculated");
        CONSOLE_LOG("Frequencies: ", signal->frequencies_axis_values[0], "..", signal->frequencies_axis_values[signal->cols-1], "; Values: ", -1.0, "..", 1.0);
        
        CONSOLE_LOG("Update signal min & max values");
        signal->amplitude_spectrum_min = 0.0;
        signal->amplitude_spectrum_max = 0.0;
        for (int i = 0; i < signal->rows; ++i)
        {
            double* amplitude_spectrum_row = signal->amplitude_spectrum + i * signal->cols;
            for (int j = 0; j < signal->cols; ++j)
            {
                if (std::abs(amplitude_spectrum_row[j]) > std::abs(signal->amplitude_spectrum_max))
                {
                    signal->amplitude_spectrum_max = amplitude_spectrum_row[j];
                }
                
                if (std::abs(amplitude_spectrum_row[j]) < std::abs(signal->amplitude_spectrum_min))
                {
                    signal->amplitude_spectrum_min = amplitude_spectrum_row[j];
                }
            }
            amplitude_spectrum_row = nullptr;
        }
        
        signal->fftw_calculated = true;
        return;
    }
	
    void CalculateThread(MSApp* app, Radar* radar, Object* object, Medium* medium, Universe* universe, SimulationRules* rules, Signal* signal)
    {
        if (!global_data.calculation_is_going)
        {
            global_data.thread_locked = true;
            global_data.calculation_is_done = false;
            global_data.calculation_is_going = true;
            global_data.todo_count = 2;
            
//            if (signal->dimensions == 2)
//            {
//                if (app->m_Image != nullptr)
//                {
//                    delete app->m_Image;
//                    app->m_Image = nullptr;
//                }
//            }
            
            global_data.todo_active_index = 1;
            global_data.percentage_status = "Сырой сигнал";
            global_data.percentage_done = 0.0;
            CalculateSignal(radar, object, medium, universe, rules, signal);
            
            if (signal->raw_calculated)
            {
                global_data.todo_active_index = 2;
                global_data.percentage_status = "Спектр сигнала";
                global_data.percentage_done = 0.0;
                CalculateSpectrum(signal, radar);
            }
            
            global_data.todo_count = 0;
            global_data.todo_active_index = 1;
            if (signal->raw_calculated && signal->fftw_calculated)
            {
                global_data.calculation_is_done = true;
                
                app->NewTask(MSAppTasks::CreateImage);
            }
            
            global_data.calculation_is_going = false;
            global_data.thread_locked = false;
        } else {
            CONSOLE_LOG("Cannot start calculation: another one is in progress");
        }
        CONSOLE_LOG("CalculateThread finished");
        return;
    }
    
    Image* CreateImage(MSApp* app, Image* image, Signal& signal, ImVec2 image_size_desired)
    {
        CONSOLE_LOG("CreateImage function start");
        printf("CreateImage image address before: %X\n", image);
        if (image != nullptr)
        {
            delete image;
            image = nullptr;
        }

        image = new Image();
        image->ready = false;
        printf("CreateImage image address now: %X\n", image);
        
        // allocate pixels memory
        uint32_t* pixels = nullptr;
        int pixels_x_max = image_size_desired.x;
        int pixels_y_max = image_size_desired.y;
        if (image_size_desired.x > signal.cols)
        {
            pixels_x_max = signal.cols;
        }
        if (image_size_desired.y > signal.rows)
        {
            pixels_y_max = signal.rows;
        }
        
        pixels = new uint32_t[pixels_x_max * pixels_y_max];
        for (size_t a = 0; a < pixels_x_max * pixels_y_max; ++a)
        {
            pixels[a] = 0xFF0000FF;
        }
        
        const int bits_red = 256*256*256;
        const int bits_green = 256*256;
        const int bits_blue = 256;
        const int bits_alpha = 1;
        int r, g, b = 0x00; // color components for a pixel
        int a = 0xFF; // alpha component for a pixel
        
        // resample data to fit current canvas size
        double signal_x_f = 0.0; // index value for scaling
        double signal_y_f = 0.0; // index value for scaling
        int signal_x = 0; // index
        int signal_y = 0; // index
        
        double signal_x_step_f = (double)signal.cols / pixels_x_max; // area to approximate
        double signal_y_step_f = (double)signal.rows / pixels_y_max; // area to approximate
        int signal_x_step = (int)round(signal_x_step_f); // area to approximate
        int signal_y_step = (int)round(signal_y_step_f); // area to approximate
        
        int signal_x_max = signal.cols;
        int signal_y_max = signal.rows;
        
        if (signal_x_step == 0)
        {
            signal_x_step_f = 1.0;
            signal_x_step = 1;
            signal_x_max = pixels_x_max;
        }
        
        if (signal_y_step == 0)
        {
            signal_y_step_f = 1.0;
            signal_y_step = 1;
            signal_y_max = pixels_y_max;
        }
        
        int signal_x_step_backup = signal_y_step;
        int signal_y_step_backup = signal_x_step;
        
        int pixels_x = 0; // index
        int pixels_y = pixels_y_max - 1; // index
        double data_temp = 0.0; // temporary storage for amplitude spectrum
        double data_max = 0.0; // what to store in current pixel
        
        CONSOLE_LOG("Imaging signal [", signal.cols, "; ", signal.rows,
                    "] with step {", signal_x_step, "; ", signal_y_step,
                    "} to (", pixels_x_max, "; ", pixels_y_max, ")");
        for (signal_x_f = 0.0; signal_x_f < (double)signal_x_max; signal_x_f += signal_x_step_f)
        {
            signal_x = (int)round(signal_x_f);
            signal_x_step = (int)round(signal_x_step_f);
            
            // make sure we don't eventually get out of array bounds
            int signal_x_remainder = signal_x_max - signal_x;
            if ((signal_x_remainder < signal_x_step) || (pixels_x == (pixels_x_max - 1)))
            {
                signal_x_step_backup = signal_x_step;
                signal_x_step = signal_x_remainder;
            }
            
            int pixels_x_index = pixels_x * pixels_y_max;
            pixels_y = pixels_y_max - 1;
            signal_y_step = signal_y_step_backup;
            
            for (signal_y_f = 0.0; signal_y_f < (double)signal_y_max; signal_y_f += signal_y_step_f)
            {
                int pixels_y_index = pixels_y * pixels_x_max;
                
                signal_y = (int)round(signal_y_f);
                signal_y_step = (int)round(signal_y_step_f);
                
                // make sure we don't eventually get out of array bounds
                int signal_y_remainder = signal_y_max - signal_y;
                if ((signal_y_remainder < signal_y_step) || (pixels_y == 0))
                {
                    signal_y_step_backup = signal_y_step;
                    signal_y_step = signal_y_remainder;
                }
                // calculate data to use for pixel color (max criteria)
                data_temp = 0.0;
                data_max = 0.0;
                for (int i = signal_x; i < (signal_x + signal_x_step); ++i)
                {
                    int i_index = i * signal.rows;
                    for (int j = signal_y; j < (signal_y + signal_y_step); ++j)
                    {
                        int j_index = j * signal.cols;
                        
                        data_temp = signal.amplitude_spectrum[j_index + i];
                        
//                        if (data_temp > data_max)
//                        {
//                            data_max = data_temp;
//                        }
                        
                        data_max += data_temp;
                    }
                }
                
                data_max = data_max / (signal_x_step * signal_y_step);
                
                // calculate pixel color from data
                app->m_ColorSchemes[app->m_ColorSchemeSelectedName].GetColor(
                            data_max / signal.amplitude_spectrum_max, 
                            app->lowerFFTWLevel, 
                            app->upperFFTWlevel, 
                            &r, &g, &b);

                *((unsigned char*)pixels + (pixels_y_index + pixels_x) * 4)     = r;
                *((unsigned char*)pixels + (pixels_y_index + pixels_x) * 4 + 1) = g;
                *((unsigned char*)pixels + (pixels_y_index + pixels_x) * 4 + 2) = b;
                *((unsigned char*)pixels + (pixels_y_index + pixels_x) * 4 + 3) = a;
                
                --pixels_y;
            }

            ++pixels_x;
        }
        CONSOLE_LOG("Pixels processed: ", pixels_x, "; ", pixels_y);
        
        // create an image
        CONSOLE_LOG("Image object creation...");
        
        if (image != nullptr)
        {
            image->Create(pixels_x_max, pixels_y_max, (void*)pixels);
            image->ready = true;
            CONSOLE_LOG("Image ready = true");
            printf("CreateImage image address: %X\n", image);
        } else {
            CONSOLE_LOG("Despite all efforts, m_Image == nullptr");
        }
        
        return image;
    }
    
	MSApp::MSApp()
	{
		CONSOLE_GREEN("Application initialized.");
		CONSOLE_LOG("Max texture size is: ", GL_MAX_TEXTURE_SIZE);
		
		SetWindowTitle("Moire Spectrum Experimentor App");
		
		m_fileDialogInfo.type = ImGuiFileDialogType_OpenFile;
		m_fileDialogInfo.title = "Открыть файл";
		m_fileDialogInfo.fileName = "";
		m_fileDialogInfo.directoryPath = std::filesystem::current_path();
		m_fileDialogInfo.fileExtensions = {
			".PNG", 
			".jpg", 
			".bmp", 
			".gif",
		};
		m_fileDialogInfo.fileExtensionSelected = -1;
		
		m_ColorSchemesMap["Greenish"] = {
			{  0,   0,   0, 255},
			{196, 255, 128, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Greenish");
		
		m_ColorSchemesMap["Matrix"] = {
			{  0,   0,   0, 255},
			{  0, 255,   0, 255},
			{  0, 255, 255, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Matrix");
		
		m_ColorSchemesMap["Black & White"] = {
			{  0,   0,   0, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Black & White");
		
		m_ColorSchemesMap["White & Black"] = {
			{255, 255, 255, 255},
			{  0,   0,   0, 255}
		};
		m_ColorSchemesNamesList.push_back("White & Black");
		
		m_ColorSchemesMap["Jet with B&W"] = {
			{  0,   0,   0, 255},
			{  0,   0, 255, 255},
			{  0, 255, 255, 255},
			{  0, 255,   0, 255},
			{255, 255,   0, 255},
			{255,   0,   0, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Jet with B&W");
		
		m_ColorSchemesMap["Jet"] = {
			{  0,   0,  64, 255},
			{  0,   0, 255, 255},
			{  0, 255, 255, 255},
			{  0, 255,   0, 255},
			{255, 255,   0, 255},
			{255,   0,   0, 255},
			{ 64,   0,   0, 255}
		};
		m_ColorSchemesNamesList.push_back("Jet");
		
		m_ColorSchemesMap["HSV"] = {
			{255,   0,   0, 255},
			{255, 255,   0, 255},
			{  0, 255,   0, 255},
			{  0, 255, 255, 255},
			{  0,   0, 255, 255},
			{255,   0, 255, 255},
			{255,   0,   0, 255}
		};
		m_ColorSchemesNamesList.push_back("HSV");
		
		m_ColorSchemesMap["Inferno"] = {
			{  0,   0,   0, 255},
			{255,   0, 255, 255},
			{255, 127, 127, 255},
			{255, 255,   0, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Inferno");
		
		m_ColorSchemesMap["Summer"] = {
			{  0, 127, 127, 255},
			{255, 255,   0, 255},
		};
		m_ColorSchemesNamesList.push_back("Summer");
		
		m_ColorSchemesMap["Seismic"] = {
			{  0,   0,   0, 255},
			{  0,   0, 255, 255},
			{255, 255, 255, 255},
			{255,   0,   0, 255},
			{127,   0,   0, 255},
		};
		m_ColorSchemesNamesList.push_back("Seismic");
		
		m_ColorSchemesMap["Hot"] = {
			{  0,   0,   0, 255},
			{255,   0,   0, 255},
			{255, 255,   0, 255},
			{255, 255, 255, 255}
		};
		m_ColorSchemesNamesList.push_back("Hot");
		
		m_ColorSchemesMap["Blue Space"] = {
			{  0,   0,   0, 255},
			{  0,   0, 128, 255},
			{  0,   0, 255, 255},
			{  0, 128, 256, 255},
			{ 64, 255,  64, 255},
			{255, 255, 255, 255},
		};
		m_ColorSchemesNamesList.push_back("Blue Space");
		
		m_ColorSchemesMap["Native"] = {
			{  0,   0,   0, 255}, // 0 // can be thrown away
			{  0,   0,  63, 255}, // 1
			{  0,   0, 127, 255},
			{  0,   0, 191, 255},
			{  0,   0, 255, 255},
			{  0,  63,   0, 255},
			{  0,  63,  63, 255},
			{  0,  63, 127, 255},
			{  0,  63, 191, 255},
			{  0,  63, 255, 255}, 
			{  0, 127,   0, 255}, // 10
			{  0, 127,  63, 255},
			{  0, 127, 127, 255},
			{  0, 127, 191, 255},
			{  0, 127, 255, 255},
			{  0, 191,   0, 255},
			{  0, 191,  63, 255},
			{  0, 191, 127, 255},
			{  0, 191, 191, 255},
			{  0, 191, 255, 255},
			{  0, 255,   0, 255}, // 20
			{  0, 255,  63, 255},
			{  0, 255, 127, 255},
			{  0, 255, 191, 255},
			{  0, 255, 255, 255},
			{ 63,   0,   0, 255},
			{ 63,   0,  63, 255},
			{ 63,   0, 127, 255},
			{ 63,   0, 191, 255},
			{ 63,   0, 255, 255},
			{ 63,  63,   0, 255}, // 30
			{ 63,  63,  63, 255},
			{ 63,  63, 127, 255},
			{ 63,  63, 191, 255},
			{ 63,  63, 255, 255},
			{ 63, 127,   0, 255},
			{ 63, 127,  63, 255},
			{ 63, 127, 127, 255},
			{ 63, 127, 191, 255},
			{ 63, 127, 255, 255},
			{ 63, 191,   0, 255}, // 40
			{ 63, 191,  63, 255},
			{ 63, 191, 127, 255},
			{ 63, 191, 191, 255},
			{ 63, 191, 255, 255},
			{ 63, 255,   0, 255},
			{ 63, 255,  63, 255},
			{ 63, 255, 127, 255},
			{ 63, 255, 191, 255},
			{ 63, 255, 255, 255},
			{127,   0,   0, 255}, // 50
			{127,   0,  63, 255},
			{127,   0, 127, 255},
			{127,   0, 191, 255},
			{127,   0, 255, 255},
			{127,  63,   0, 255},
			{127,  63,  63, 255},
			{127,  63, 127, 255},
			{127,  63, 191, 255},
			{127,  63, 255, 255},
			{127, 127,   0, 255}, // 60
			{127, 127,  63, 255},
			{127, 127, 127, 255},
			{127, 127, 191, 255},
			{127, 127, 255, 255},
			{127, 191,   0, 255},
			{127, 191,  63, 255},
			{127, 191, 127, 255},
			{127, 191, 191, 255},
			{127, 191, 255, 255},
			{127, 255,   0, 255}, // 70
			{127, 255,  63, 255},
			{127, 255, 127, 255},
			{127, 255, 191, 255},
			{127, 255, 255, 255},
			{191,   0,   0, 255},
			{191,   0,  63, 255},
			{191,   0, 127, 255},
			{191,   0, 191, 255},
			{191,   0, 255, 255},
			{191,  63,   0, 255}, // 80
			{191,  63,  63, 255},
			{191,  63, 127, 255},
			{191,  63, 191, 255},
			{191,  63, 255, 255},
			{191, 127,   0, 255},
			{191, 127,  63, 255},
			{191, 127, 127, 255},
			{191, 127, 191, 255},
			{191, 127, 255, 255},
			{191, 191,   0, 255}, // 90
			{191, 191,  63, 255},
			{191, 191, 127, 255},
			{191, 191, 191, 255},
			{191, 191, 255, 255},
			{191, 255,   0, 255},
			{191, 255,  63, 255},
			{191, 255, 127, 255},
			{191, 255, 191, 255},
			{191, 255, 255, 255},
			{255,   0,   0, 255}, // 100
			{255,   0,  63, 255},
			{255,   0, 127, 255},
			{255,   0, 191, 255},
			{255,   0, 255, 255},
			{255,  63,   0, 255},
			{255,  63,  63, 255},
			{255,  63, 127, 255},
			{255,  63, 191, 255},
			{255,  63, 255, 255},
			{255, 127,   0, 255}, // 110
			{255, 127,  63, 255},
			{255, 127, 127, 255},
			{255, 127, 191, 255},
			{255, 127, 255, 255},
			{255, 191,   0, 255},
			{255, 191,  63, 255},
			{255, 191, 127, 255},
			{255, 191, 191, 255},
			{255, 191, 255, 255},
			{255, 255,   0, 255}, // 120
			{255, 255,  63, 255},
			{255, 255, 127, 255},
			{255, 255, 191, 255},
			{255, 255, 255, 255}, // 124
		};
		m_ColorSchemesNamesList.push_back("Native");
		
		m_ColorSchemeSelected = 0;
		m_ColorSchemeSelectedName = m_ColorSchemesNamesList[m_ColorSchemeSelected];
		
		for (auto it = m_ColorSchemesMap.begin(); it != m_ColorSchemesMap.end(); ++it)
		{
			m_ColorSchemes[it->first] = {it->first, {}};
			for (int i = 0; i < m_ColorSchemesMap[it->first].size(); ++i)
			{
				m_ColorSchemes[it->first].colors.push_back(m_ColorSchemesMap[it->first][i]);
			}
			CONSOLE_LOG(it->first, " size: ", m_ColorSchemes[it->first].colors.size());
		}
		
		m_CurrentMode = MSAppMode::Idle;
		CONSOLE_LOG("Enter Idle mode");
		
		NewTask(MSAppTasks::Idle);
		CONSOLE_LOG("Add a new Task: Idle");
		
		// setup File Dialog initial values
		m_fileDialogInfo.type = ImGuiFileDialogType_OpenFile;
		m_fileDialogInfo.title = "Открыть файл";
		m_fileDialogInfo.fileName = "";
		m_fileDialogInfo.directoryPath = std::filesystem::current_path();
		m_fileDialogInfo.fileExtensions = {
			".PNG", 
			".jpg", 
			".bmp", 
			".gif",
		};
		m_fileDialogInfo.fileExtensionSelected = -1;
        
        // setup radar
        radar.transmit_power = 1.0;
        radar.starting_frequency = 10e9;
        radar.modulation_time = 0.001;
        radar.deviation_frequency = 30e6;
        radar.sampling_frequency = radar.deviation_frequency * 2;
        radar.time_sample_width = 1/radar.sampling_frequency;
        radar.observation_time = 0.001;
        radar.position.x = 0.0;
        radar.position.y = 0.0;
        radar.position.z = 1000.0;
        
        // setup object
        object.position.x = 0.0;
        object.position.y = 1000.0;
        object.position.z = 0.0;
        object.linear_size = 50000.0;
        object.distance_between_elements = 1.0;
        object.size = (int)(round(object.linear_size / object.distance_between_elements));
        object.structure.resize(object.size);
        for (int i = 0; i < object.size; ++i)
        {
            object.structure[i].esr = 1.0;
            object.structure[i].phase_shift = 0.0;
            object.structure[i].position.x = 0.0;
            object.structure[i].position.y = 0.0 + i * object.distance_between_elements;
            object.structure[i].position.z = 0.0;
//            object.structure[i].position.z = 0.0 + ((double)rand() / RAND_MAX - 0.5);
        }
        CONSOLE_LOG("Object consists of ", object.size, " elements");
        
        // setup universe
        // actually, default values go well
        
        // setup medium
        medium.propagation_coefficient = 1.0;
        
        // setup signal
//        signal.cols = (int)(round(radar.observation_time / radar.time_sample_width));
//        signal.rows = 1;
//        signal.dimensions = 1;
//        global_data.signal_dimensions_mode = SignalDimensionsMode::Time_Amplitude;
        
        // setup simulation rules
        rules.no_amplitude = true;
        rules.one_distance = true;
        rules.stop_and_go = true;
        rules.spare_memory_on_frequencies = true;
        
	}
	
	MSApp::~MSApp()
	{
		delete m_Image;
		m_Image = nullptr;
        
//        delete[] signal.raw_data;
//        signal.raw_data = nullptr;
		
		CONSOLE_GREEN("Application is shut down.");
	}
	
	void MSApp::SetupResources()
	{
		// image stuff
//		m_Image = new Image();
//		m_Image->Load("data/img/temp.png");
	}
	
	void MSApp::Logic()
	{
		while (m_TaskStack.size() > 0)
		{
			m_CurrentTask = m_TaskStack.back();
			m_TaskStack.pop_back();
			
			switch (m_CurrentTask)
			{
    			case MSAppTasks::Idle:
				// do nothing
    				break;
    			case MSAppTasks::CalculateFFT:
    				{
                        if (signal.raw_calculated)
                        {
                            signal.Reset();
                        }
                        
                        std::thread calculation(CalculateThread, this, &radar, &object, &medium, &universe, &rules, &signal);
                        calculation.detach();
                        CONSOLE_LOG("CalculateFFT task complete");
    				}
    				break;
                case MSAppTasks::CreateImage:
                    {
                        printf("m_Image address: %X\n", m_Image);
                        CONSOLE_LOG("CreateImage start");
                        m_Image = CreateImage(this, m_Image, signal, {m_Image_width, m_Image_height});
                        CONSOLE_LOG("Image created");
                        printf("m_Image address: %X\n", m_Image);
                        if (m_Image != nullptr)
                        {
                            printf("m_Image->ready: %X\n", m_Image->ready);
                        }
                    }
                    break;
    			case MSAppTasks::LoadImage:
    				{
    					CONSOLE_LOG("Loading image from path: ", m_fileDialogInfo.resultPath.string());
    					m_Image->Unload();
    					m_Image->Load(m_fileDialogInfo.resultPath.string());
    					
    					lowerFFTWLevel = m_Image->fftw_data_min;
    					
    					m_FFTCalculated = false;
    				}
    				break;
    			case MSAppTasks::SaveFFT:
    				{
    					CONSOLE_LOG("Saving FFT image to path: ", m_fileDialogInfo.resultPath.string());
    					m_Image->Save(m_fileDialogInfo.resultPath.string());
    				}
    				break;
    			case MSAppTasks::Exit:
    				{
    					doExit = true;
    				}
    				break;
    			default:
    				break;
			}
		}
		m_CurrentTask = MSAppTasks::Idle;
	}
	
	void MSApp::GUIContent() 
	{
		bool p_open = true;
		static bool use_work_area = true;
		// static ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
		static ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoSavedSettings;
		
		TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x * m_WindowScale.x;
		TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing() * m_WindowScale.y;
//		CONSOLE_LOG("Window scale: (", m_WindowScale.x, ", ", m_WindowScale.y, ")");
		
		// ====================================================================================
		// The Application starts here
		// ------------------------------------------------------------------------------------
		
		ShowMainMenu();
		
		// We demonstrate using the full viewport area or the work area (without menu-bars, task-bars etc.)
		// Based on your use case you may want one or the other.
		ImGuiViewport* viewport = ImGui::GetMainViewport();
		ImGui::SetNextWindowPos(use_work_area ? viewport->WorkPos : viewport->Pos);
		ImGui::SetNextWindowSize(use_work_area ? viewport->WorkSize: viewport->Size);
		// ------------------------------------------------------------------------------------
		
		if (ImGui::Begin("Example: Fullscreen window", &p_open, flags))
		{
            if (ImGui::BeginTabBar("TrackSignalModelTabBar", ImGuiTabBarFlags_None))
            {
                if (ImGui::BeginTabItem("Демонстрация сигнала##tabitem"))
                {
                    if (ImGui::BeginTable("Демонстрация сигнала##table", 2, 0, {ImGui::GetContentRegionAvail().x, 0.0}))
                    {
                        static float columnWidth1 = 3 * m_WindowCurrentWidth / 4; // 960
                        static float columnWidth2 = 1 * m_WindowCurrentWidth / 4; // 320
                        
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, columnWidth1 * m_WindowScale.x);
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, columnWidth2 * m_WindowScale.x);
                        
                        ImGui::TableNextColumn();
                        static float columnHeight = ImGui::GetContentRegionAvail().y - 20;

                        if ((signal.dimensions == 1) && signal.fftw_calculated)
                        {
                            if (ImPlot::BeginPlot("Амплитудно-частотный спектр преобразованного сигнала", ImGui::GetContentRegionAvail()))
                            {
                                if (global_data.calculation_is_done)
                                {
                                    switch (m_FourierSpectrumMode)
                                    {
                                        case FourierSpectrumMode::Amplitude:
                                            {
                                                ImPlot::SetupAxes("Частота биений, Гц", "Вт * с", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_NoTickLabels);
        //                                        ImPlot::SetupAxisLimits(ImAxis_X1, signal.frequencies_axis_values[0], signal.frequencies_axis_values[signal.cols-1]);
                                                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Linear);
                                                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Linear);
        //                                        ImPlot::SetupAxisTicks(ImAxis_Y1, 0.0, 8000, 8, nullptr, false);
                                                ImPlot::PushStyleColor(ImPlotCol_Line, {1.0, 0.6, 0.1, 1.0});
                                                ImPlot::PlotLine("Радиолокационный профиль", signal.frequencies_axis_values, signal.amplitude_spectrum, signal.cols);
                                            }
                                            break;
                                        case FourierSpectrumMode::Phase:
                                            {
                                                ImPlot::SetupAxes("Частота биений, Гц", "1", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_NoTickLabels);
    //                                            ImPlot::SetupAxisLimits(ImAxis_X1, signal.frequencies_axis_values[0], signal.frequencies_axis_values[signal.cols-1]);
                                                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Linear);
                                                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Linear);
    //                                            ImPlot::SetupAxisTicks(ImAxis_Y1, -pi, pi, 8);
                                                ImPlot::PushStyleColor(ImPlotCol_Line, {1.0, 0.6, 0.1, 1.0});
                                                ImPlot::PlotLine("Радиолокационный профиль", signal.frequencies_axis_values, signal.phase_spectrum, signal.cols);
                                            }
                                            break;
                                    }
                                } else {
//                                    ImPlot::PlotDummy("###JustAPlot");
                                }
                                ImPlot::EndPlot();
                            }
                        }
                        
                        if (signal.dimensions == 2)
                        {
                            // draw a 2d-plot which is a heatmap
                            if (m_Image != nullptr)
                            {
                                if (m_Image->ready)
                                {
                                    ImGui::Image((ImTextureID)(m_Image->texture->GetID()),
                                                 ImVec2(columnWidth1*m_WindowScale.x, 
                                                        columnHeight*m_WindowScale.y),
                                                 ImVec2(0, 0),
                                                 ImVec2(1, 1));
                            
//                                    // DEBUG SECTION
//                                    if (ImPlot::BeginPlot("Амплитудно-частотный спектр преобразованного сигнала", ImGui::GetContentRegionAvail()))
//                                    {
//                                        ImPlot::SetupAxes("Частота биений, Гц", "Вт * с", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_NoTickLabels);
//                                        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Linear);
//                                        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Linear);
//                                        ImPlot::PushStyleColor(ImPlotCol_Line, {1.0, 0.6, 0.1, 1.0});
//                                        ImPlot::PlotLine("Радиолокационный профиль", signal.frequencies_axis_values, signal.amplitude_spectrum + 999 * signal.cols, signal.cols/2);
//                                        ImPlot::EndPlot();
//                                    }
                                }
                            }
                        }
                        
                        ImGui::TableNextColumn();
                        
                        ImGui::SeparatorText("Настройка усиления и цветовой схемы");
                        ImGui::Text("Нижний порог: ");
                        ImGui::SameLine(150);
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputFloat("###lowerFFT", &lowerFFTWLevel, 0.01, 0.1, "%.10f"))
                        {
                            if (lowerFFTWLevel < 0.0)
                            {
                                lowerFFTWLevel = 0.0;
                            }
                            if (lowerFFTWLevel > 1.0)
                            {
                                lowerFFTWLevel = 1.0;
                            }
                            NewTask(MSAppTasks::CreateImage);
                        }
                        
                        ImGui::Text("Верхний порог: ");
                        ImGui::SameLine(150);
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputFloat("###upperFFT", &upperFFTWlevel, 0.01, 0.1, "%.10f"))
                        {
                            if (upperFFTWlevel < 0.0)
                            {
                                upperFFTWlevel = 0.0;
                            }
                            if (upperFFTWlevel > 1.0)
                            {
                                upperFFTWlevel = 1.0;
                            }
                            NewTask(MSAppTasks::CreateImage);
                        }
                        
                        ImGui::Text("Раскраска: ");
                        ImGui::SameLine(150);
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::BeginCombo("###ColorSchemeSelector", m_ColorSchemesNamesList[m_ColorSchemeSelected].c_str()))
                        {
                            static bool isSelected = false;
                            for (int i = 0; i < m_ColorSchemesNamesList.size(); ++i)
                            {
                                isSelected = (i == m_ColorSchemeSelected);
                                if (ImGui::Selectable(m_ColorSchemesNamesList[i].c_str(), isSelected))
                                {
                                    m_ColorSchemeSelected = i;
                                    m_ColorSchemeSelectedName = m_ColorSchemesNamesList[i];
                                    NewTask(MSAppTasks::CreateImage);
                                }
                            }
                            ImGui::EndCombo();
                        }
                        
                        ImGui::Separator();
                        ImGui::Text("Вид частотного спектра:"); 
                        ImGui::SameLine();
                        switch (m_FourierSpectrumMode)
                        {
                            case FourierSpectrumMode::Amplitude:
                                {
                                    ImGui::TextColored({1.0, 0.8, 0.2, 1.0}, "Амплитудный");
                                }
                                break;
                            case FourierSpectrumMode::Phase:
                                {
                                    ImGui::TextColored({1.0, 0.8, 0.2, 1.0}, "Фазовый");
                                }
                                break;
                        }
                        
                        ImGui::PushStyleColor(ImGuiCol_Button, {0.5, 0, 0.25, 1.0});
                        if (ImGui::Button("Амплитудный", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            m_FourierSpectrumMode = FourierSpectrumMode::Amplitude;
                            if (!signal.fftw_calculated)
                            {
                                CONSOLE_LOG("new task: CalculateFFT");
                                NewTask(MSAppTasks::CalculateFFT);
                            } else {
//					m_Image->UpdatePixels(&m_ColorSchemes[m_ColorSchemeSelectedName], m_FourierSpectrumMode);
//					m_Image->UpdateTexture();
                            }
                        }
                        ImGui::PopStyleColor();
                        
                        ImGui::SameLine(170 * m_WindowScale.x);
                        
                        ImGui::PushStyleColor(ImGuiCol_Button, {0.25, 0, 0.5, 1.0});
                        if (ImGui::Button("Фазовый", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            m_FourierSpectrumMode = FourierSpectrumMode::Phase;
                            if (!signal.fftw_calculated)
                            {
                                CONSOLE_LOG("new task: CalculateFFT");
                                NewTask(MSAppTasks::CalculateFFT);
                            } else {
//					m_Image->UpdatePixels(&m_ColorSchemes[m_ColorSchemeSelectedName], m_FourierSpectrumMode);
//					m_Image->UpdateTexture();
                            }
                        }
                        ImGui::PopStyleColor();
                        
                        ImGui::PushStyleColor(ImGuiCol_Button, {0.25, 0, 0.25, 1.0});
                        if (ImGui::Button("Пересчитать", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            CONSOLE_LOG("new task: CalculateFFT");
                            NewTask(MSAppTasks::CalculateFFT);
                        }
                        ImGui::PopStyleColor();
                        
                        ImGui::SameLine(170 * m_WindowScale.x);
                        
                        ImGui::PushStyleColor(ImGuiCol_Button, {0.25, 0.25, 0.5, 1.0});
                        if (ImGui::Button("Просто кнопка", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                        }
                        ImGui::PopStyleColor();
                        
                        ImGui::SeparatorText("Режим симуляции");
                        if (ImGui::Button("U(t) | P(fб)", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            global_data.signal_dimensions_mode = SignalDimensionsMode::Time_Amplitude;
                            signal.dimensions = 1;
                            
                            CONSOLE_LOG("new task: CalculateFFT");
                            NewTask(MSAppTasks::CalculateFFT);
                        }
                        
                        if (ImGui::Button("S(f0, fб)", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            global_data.signal_dimensions_mode = SignalDimensionsMode::BeatFrequency_CarrierFrequency;
                            global_data.signal_rows_variable_min = 0e6;
                            global_data.signal_rows_variable_max = 10e9;
                            global_data.signal_rows_variable_step = 1e7;
                            signal.dimensions = 2;
                            
                            CONSOLE_LOG("new task: CalculateFFT");
                            NewTask(MSAppTasks::CalculateFFT);
                        }
                        
                        ImGui::SameLine(170 * m_WindowScale.x);
                        if (ImGui::Button("S(f0, fб), var(B)", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            global_data.signal_dimensions_mode = SignalDimensionsMode::BeatFrequency_CarrierFrequency_B;
                            signal.dimensions = 2;
                            
                            CONSOLE_LOG("new task: CalculateFFT");
                            NewTask(MSAppTasks::CalculateFFT);
                        }
                        
                        if (ImGui::Button("S(h, fб)", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            global_data.signal_dimensions_mode = SignalDimensionsMode::BeatFrequency_Height;
                            signal.dimensions = 2;
                            global_data.signal_rows_variable_min = 0;
                            global_data.signal_rows_variable_max = 2000;
                            global_data.signal_rows_variable_step = 2e0;
                            
                            CONSOLE_LOG("new task: CalculateFFT");
                            NewTask(MSAppTasks::CalculateFFT);
                        }
                        
                        if (ImGui::Button("S(dy, fб)", {150 * m_WindowScale.x, 30 * m_WindowScale.y}))
                        {
                            global_data.signal_dimensions_mode = SignalDimensionsMode::BeatFrequency_Dy;
                            signal.dimensions = 2;
                            if (!signal.fftw_calculated)
                            {
                                CONSOLE_LOG("new task: CalculateFFT");
                                NewTask(MSAppTasks::CalculateFFT);
                            } else {
//					m_Image->UpdatePixels(&m_ColorSchemes[m_ColorSchemeSelectedName], m_FourierSpectrumMode);
//					m_Image->UpdateTexture();
                            }
                        }
                        
                        ImGui::SeparatorText("Статус");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (global_data.todo_count > 0)
                        {
                            ImGui::Text("(%d/%d) - %s - %.1f%%", 
                                        global_data.todo_active_index, 
                                        global_data.todo_count, 
                                        global_data.percentage_status.c_str(), 
                                        global_data.percentage_done);
                        } else {
                            ImGui::Text("Всё в порядке");
                        }
                        
                        ImGui::EndTable();
                    }
                    
                    ImGui::EndTabItem();
                }
                
                if (ImGui::BeginTabItem("Настройка симуляции"))
                {
                    static float columnWidth2 = 1 * m_WindowCurrentWidth / 4; // 320
                    static float columnHeight = ImGui::GetContentRegionAvail().y - 20;
                    
                    ImGui::SeparatorText("Настройка симуляции##192810293812");
                    if (ImGui::BeginTable("Радар и сигнал", 2, 0, {ImGui::GetContentRegionAvail().x, 0.0}))
                    {
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, ImGui::GetContentRegionAvail().x/2);
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, ImGui::GetContentRegionAvail().x/2);
                        
                        ImGui::TableNextColumn();
                        ImGui::SeparatorText("Радар и сигнал");
                        ImGui::Text("Начальная частота, Гц");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarStartingFrequency", &radar.starting_frequency, 10000, 100000, "%.1f"))
                        {
                        }
                        
                        ImGui::Text("Время наблюдения, с");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarObservationTime", &radar.observation_time, 0.000001, 0.001, "%.6f"))
                        {
                        }
                        
                        ImGui::Text("Интервал модуляции, с");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarModulationTime", &radar.modulation_time, 0.000001, 0.001, "%.6f"))
                        {
                        }
                        
                        ImGui::Text("Девиация частоты, Гц/с");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarDeviationFrequency", &radar.deviation_frequency, 10000, 100000, "%.1f"))
                        {
                            if (rules.spare_memory_on_frequencies)
                            {
                                rules.beat_frequency_min = 0.0;
                                
                                // this memory saving technique is based on the assumption that the object is basically a line perpendicular to the X axis
                                double distance_max = sqrt(pow(radar.position.x - (object.position.x + object.structure[object.size - 1].position.x), 2) +
                                                           pow(radar.position.y - (object.position.y + object.structure[object.size - 1].position.y), 2) +
                                                           pow(radar.position.z - (object.position.z + object.structure[object.size - 1].position.z), 2));
                                double tk_max = 0.5*distance_max / universe.c;
                                rules.beat_frequency_max = radar.deviation_frequency * tk_max;
                                radar.sampling_frequency = 2 * rules.beat_frequency_max;
                                radar.time_sample_width = 1/radar.sampling_frequency;
                                CONSOLE_LOG("Sparing memory mode ON");
                                CONSOLE_LOG("tk_max: ", tk_max, " sec");
                                CONSOLE_LOG("beat frequency max: ", rules.beat_frequency_max, " Hz");
                                CONSOLE_LOG("sampling frequency: ", radar.sampling_frequency, " Hz");
                            } else {
                                // still spare memory
                                rules.beat_frequency_max = radar.deviation_frequency * radar.observation_time;
                                radar.sampling_frequency = 2 * rules.beat_frequency_max;
                                radar.time_sample_width = 1/radar.sampling_frequency;
                                CONSOLE_LOG("Sparing memory mode OFF");
                                CONSOLE_LOG("beat frequency max: ", rules.beat_frequency_max, " Hz");
                                CONSOLE_LOG("sampling frequency: ", radar.sampling_frequency, " Hz");
                            }
                        }
                        
                        ImGui::Text("Частота дискретизации, Гц");
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarSamplingFrequency", &radar.sampling_frequency, 10000, 100000, "%.1f"))
                        {
                            radar.time_sample_width = 1/radar.sampling_frequency;
                        }
                        
                        ImGui::Text("Расположение радара, м");
                        ImGui::Text("X: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarPositionX", &radar.position.x, 0.1, 1, "%.1f"))
                        {
                        }
                        
                        ImGui::Text("Y: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarPositionY", &radar.position.y, 0.1, 1, "%.1f"))
                        {
                        }
                        
                        ImGui::Text("Z: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###RadarPositionZ", &radar.position.z, 0.1, 1, "%.1f"))
                        {
                        }
                        
                        ImGui::TableNextColumn();
                        ImGui::SeparatorText("Объект");
                        ImGui::Text("Линейный размер, м");
                        if (ImGui::InputDouble("###ObjectLinearSize", &object.linear_size, 1.0, 10.0, "%.2f"))
                        {
                            object.size = (int)(round(object.linear_size / object.distance_between_elements));
                            object.structure.resize(object.size);
                            for (int i = 0; i < object.size; ++i)
                            {
                                object.structure[i].esr = 1.0;
                                object.structure[i].phase_shift = 0.0;
                                object.structure[i].position.x = 0.0;
                                object.structure[i].position.y = 0.0 + i * object.distance_between_elements;
                                object.structure[i].position.z = 0.0;
                            }
                            CONSOLE_LOG("Object consists of ", object.size, " elements");
                        }
                        ImGui::Text("Расстояние между элементами, м");
                        if (ImGui::InputDouble("###ObjectElementLinearSize", &object.distance_between_elements, 1.0, 10.0, "%.2f"))
                        {
                            if (object.distance_between_elements != 0.0)
                            {
                                object.size = (int)(round(object.linear_size / object.distance_between_elements));
                                object.structure.resize(object.size);
                                for (int i = 0; i < object.size; ++i)
                                {
                                    object.structure[i].esr = 1.0;
                                    object.structure[i].phase_shift = 0.0;
                                    object.structure[i].position.x = 0.0;
                                    object.structure[i].position.y = 0.0 + i * object.distance_between_elements;
                                    object.structure[i].position.z = 0.0;
                                }
                                CONSOLE_LOG("Object consists of ", object.size, " elements");
                            }
                        }
                        ImGui::Text("Расположение объекта, м");
                        ImGui::Text("X: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###ObjectPositionX", &object.position.x, 0.1, 1, "%.1f"))
                        {
                        }
                        
                        ImGui::Text("Y: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###ObjectPositionY", &object.position.y, 0.1, 1, "%.1f"))
                        {
                        }
                        
                        ImGui::Text("Z: ");
                        ImGui::SameLine();
                        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
                        if (ImGui::InputDouble("###ObjectPositionZ", &object.position.z, 0.1, 1, "%.1f"))
                        {
                        }
                        ImGui::EndTable();
                    }
                    
                    
                    if (ImGui::BeginTable("Настройки симуляции", 2, 0, {ImGui::GetContentRegionAvail().x, 0.0}))
                    {
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, 20*m_WindowScale.x);
                        ImGui::TableSetupColumn("###", ImGuiTableColumnFlags_WidthFixed, (columnWidth2 - 20)*m_WindowScale.x);
                        ImGui::TableNextColumn();
                        if (ImGui::Checkbox("###NoAmplitude", &rules.no_amplitude))
                        {
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("Без учёта амплитуд");
                        
                        ImGui::TableNextColumn();
                        if (ImGui::Checkbox("###OneDistance", &rules.one_distance))
                        {
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("Постоянная дальность");
                        
                        ImGui::TableNextColumn();
                        if (ImGui::Checkbox("###SpareMemory", &rules.spare_memory_on_frequencies))
                        {
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("Экономить память");
                        
                        ImGui::TableNextColumn();
                        if (ImGui::Checkbox("###StopAndGo", &rules.stop_and_go))
                        {
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("Stop & Go");
                        ImGui::EndTable();
                    }
                    
                    ImGui::EndTabItem();
                }
                
                ImGui::EndTabBar();
            }
			ImGui::End();
		}
		// ====================================================================================
	}
	
	void MSApp::NewTask(MSAppTasks task)
	{
		m_ErrorFlag = false;
		m_ErrorMessage = "";
		m_TaskStack.push_back(task);
	}
	
	void MSApp::ShowMainMenu()
	{
		if (ImGui::BeginMainMenuBar())
		{
			if (ImGui::BeginMenu("Меню"))
			{
				if (ImGui::MenuItem("Открыть изображение...")) 
				{
					
					CONSOLE_LOG("Ask to open popup");
					m_fileDialogInfo.type = ImGuiFileDialogType_OpenFile;
					m_fileDialogInfo.title = "Открыть файл";
					m_fileDialogInfo.fileName = "";
					m_fileDialogInfo.fileExtensions = {
						".PNG", 
						".jpg", 
						".bmp", 
						".gif",
					};
					m_fileDialogInfo.fileExtensionSelected = -1;
					m_fileDialogOpen = true;
				}

				ImGui::Separator();
				if (ImGui::MenuItem("Сохранить изображение..."))
				{
					CONSOLE_LOG("Ask to open popup");
					m_fileDialogInfo.type = ImGuiFileDialogType_SaveFile;
					m_fileDialogInfo.title = "Сохранить файл";
					m_fileDialogInfo.fileName = "fft_image.png";
					m_fileDialogInfo.fileExtensions = {
						".png", 
						".jpg", 
						".bmp", 
						".gif",
					};
					m_fileDialogInfo.fileExtensionSelected = -1;
					m_fileDialogOpen = true;
					CONSOLE_LOG("Add a new Task: Save FFT image");
				}
				ImGui::Separator();
				if (ImGui::MenuItem("Выход")) 
				{
					NewTask(MSAppTasks::Exit);
					CONSOLE_LOG("Add a new Task: Exit");
				}
				ImGui::EndMenu();
			}
			
			if (ImGui::MenuItem("Фурье-образ"))
			{
				NewTask(MSAppTasks::CalculateFFT);
				CONSOLE_LOG("Add a new Task: Calculate Fourier Transform");
			}
			
			ImGui::EndMainMenuBar();
		}
		
		if (m_fileDialogOpen)
		{
			switch (m_fileDialogInfo.type)
			{
			case ImGuiFileDialogType_OpenFile:
				{
					if (ImGui::FileDialog(&m_fileDialogOpen, &m_fileDialogInfo))
					{
						// Result path in: m_fileDialogInfo.resultPath
						NewTask(MSAppTasks::LoadImage);
						CONSOLE_LOG("Add a new Task: Load an image");
						m_fileDialogOpen = false;
						ImGui::CloseCurrentPopup();
					}
				}
				break;
			case ImGuiFileDialogType_SaveFile:
				{
					if (ImGui::FileDialog(&m_fileDialogOpen, &m_fileDialogInfo))
					{
						// Result path in: m_fileDialogInfo.resultPath
						NewTask(MSAppTasks::SaveFFT);
						CONSOLE_LOG("Add a new Task: Save an FFT image");
						m_fileDialogOpen = false;
						ImGui::CloseCurrentPopup();
					}
				}
				break;
			}
		}
	}
	
	void MSApp::ShowContent()
	{
	}
	
	App* CreateApplication()
	{
		MSApp* app = new MSApp();
		app->SetFPS(SAVANNAH_FPS60);
		return app;
	}
}
