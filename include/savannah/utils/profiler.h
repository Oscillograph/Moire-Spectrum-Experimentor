#ifndef SAVANNAH_PROFILER_H
#define SAVANNAH_PROFILER_H

#include <savannah/proto-core.h>
#include <savannah/utils/time.h>
#include <thread>

// Most of this Profiler's code has been taken from Instrumentor.h file from:
// https://gist.github.com/TheCherno/31f135eea6ee729ab5f26a6908eb3a5e/raw/55619bb8dd4349895f7d5f2520cc35a0e6d38c09/Instrumentor.h

#ifdef SAVANNAH_PROFILING
	// Profiling to files in json format. These macros' order represents usage in a project.
	// - each session is saved in a different file;
	// - "Scope" means a timer shall auto-destroy itself once its scope ends, so it might be a good practice to place the corresponding code in brackets;
	// - "Function" is an advanced version of "Scope" timers, it should be places at the start of the function
	#define SAVANNAH_PROFILE_SESSION_START(name, filepath) 	::BSE::Profiler::BeginSession(name, filepath);
	#define SAVANNAH_PROFILE_SCOPE(name)					::BSE::ProfilerTimer timer##__LINE__(name);
	#define SAVANNAH_PROFILE_FUNCTION() 					::BSE::ProfilerTimer funcTimer##__LINE__(__PRETTY_FUNCTION__);
	#define SAVANNAH_PROFILE_SESSION_END() 					::BSE::Profiler::EndSession();
			
	// A simpler version: moment profiling which can be read by ImGuiLayers, etc.
	#define SAVANNAH_CORE_PROFILER(name)	Timer timer##__LINE__(name, [&](BSE::TimerResult result){ BSE::Profiler::Push(result); });
	#define SAVANNAH_PROFILER(name)			Timer timer##__LINE__(name, [&](BSE::TimerResult result){ BSE::Profiler::Push(result); });
#else
	// profiling to files in json format
	#define SAVANNAH_PROFILE_SESSION_START(name, filepath)
	#define SAVANNAH_PROFILE_SCOPE(name)				
	#define SAVANNAH_PROFILE_FUNCTION()					
	#define SAVANNAH_PROFILE_SESSION_END() 				
			
	// a simpler version: moment profiling
	#define SAVANNAH_CORE_PROFILER(name)	
	#define SAVANNAH_PROFILER(name)	
#endif

namespace Savannah {
	struct ProfileResult {
		std::string Name;
		long long Start; 
		long long End;
		uint32_t ThreadID;
	};
	
	struct ProfilerSession {
		std::string Name;
	};
	
	// Usage example of a "simpler profiler".
	// 
	// SOMEWHERE:
	//	SAVANNAH_PROFILER(u8"my_operation::step");
	//
	// SOMEWHERE-IN-GUI-MANAGEMENT:
	// 	auto profileResults = Profiler::Read();
	// 	for (auto kv : profileResults){
	//		std::string name = kv.first;
	// 		float time = kv.second;
	//		MY_SHOW_RESULTS_FUNCTION(name, time);
	// 	}
	
	class SAVANNAH_API Profiler {
	public:
		inline static void Push(TimerResult result){
			m_ProfileResults[result.Name] = result.Time;
		}
		inline static std::unordered_map<std::string, float> Read(){
			return m_ProfileResults;
		}
		inline static void Flush(){
			m_ProfileResults.clear();
		}
		inline static ProfilerSession* GetSession(){ 
			return m_CurrentSession; 
		}
		inline static void BeginSession(std::string name, std::string filepath = "./debug/runtime.json"){
			m_OutputStream.open(filepath.c_str());
			m_Filepath = filepath;
			WriteHeader();
			m_CurrentSession = new ProfilerSession{name};
		}
		
		inline static void EndSession() {
			WriteFooter();
			
			m_OutputStream.close();
			delete m_CurrentSession;
			m_CurrentSession = nullptr;
			m_ProfileCount = 0;
		}
		
		inline static void WriteHeader(){
			m_OutputStream << "{\"otherData\": {},\"traceEvents\":[";
			m_OutputStream.flush();
		}
		
		inline static void WriteProfile(const ProfileResult& result)
		{
			if (m_ProfileCount++ > 0)
				m_OutputStream << ",";
			
			std::string name = result.Name;
			std::replace(name.begin(), name.end(), '"', '\'');
			
			m_OutputStream << "{";
			m_OutputStream << "\"cat\":\"function\",";
			m_OutputStream << "\"dur\":" << (result.End - result.Start) << ',';
			m_OutputStream << "\"name\":\"" << name << "\",";
			m_OutputStream << "\"ph\":\"X\",";
			m_OutputStream << "\"pid\":0,";
			m_OutputStream << "\"tid\":" << result.ThreadID << ",";
			m_OutputStream << "\"ts\":" << result.Start;
			m_OutputStream << "}";
			
			m_OutputStream.flush();
		}
		
		inline static void WriteFooter() {
			m_OutputStream << "]}";
			m_OutputStream.flush();
		}
		
	private:
		static std::unordered_map<std::string, float> m_ProfileResults;
		static ProfilerSession* m_CurrentSession;
		static std::string m_Filepath;
		static std::ofstream m_OutputStream;
		static int m_ProfileCount;
	}; 
	
	// Special Profiler's version of timer
	class SAVANNAH_API ProfilerTimer {
	public:
		ProfilerTimer(std::string name)
		: m_Name(name), m_Stopped(false)
		{
			m_StartTimePoint = std::chrono::high_resolution_clock::now();
		}
		
		~ProfilerTimer(){
			if (!m_Stopped){
				Stop();
			}
		}
		
		void Stop(){
			auto endTimePoint = std::chrono::high_resolution_clock::now();
			long long start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimePoint).time_since_epoch().count();
			long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimePoint).time_since_epoch().count();
			m_Stopped = true;
			
			// float duration = (end - start) * 0.001f;
			
			uint32_t threadID = std::hash<std::thread::id>{}(std::this_thread::get_id());
			Profiler::WriteProfile({ m_Name, start, end, threadID });
		}
		
	private:
		std::string m_Name;
		std::chrono::high_resolution_clock::time_point m_StartTimePoint;
		bool m_Stopped;
	};
}
#endif
