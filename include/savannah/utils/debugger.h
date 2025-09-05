#ifndef SAVANNAH_DEBUGGER_H
#define SAVANNAH_DEBUGGER_H

#include <savannah/proto-core.h>

#ifdef SAVANNAH_SPDLOG
	#include <savannah/external/spdlog/spdlog.h>
	#include <savannah/external/spdlog/sinks/stdout_color_sinks.h>
	#include <savannah/external/spdlog/fmt/ostr.h>
	
	namespace Savannah{
		class SAVANNAH_API Debugger {
		public:
			static void Init();
			inline static std::shared_ptr<spdlog::logger>& GetCoreDebugger(){
				return s_CoreDebugger;
			};
			inline static std::shared_ptr<spdlog::logger>& GetClientDebugger(){
				return s_ClientDebugger;
			};
		private:
			static std::shared_ptr<spdlog::logger> s_CoreDebugger;
			static std::shared_ptr<spdlog::logger> s_ClientDebugger;
		};
	}
	
	// Savannah core log macros
	#define SAVANNAH_CORE_TRACE(...)		::Savannah::Debugger::GetCoreDebugger()->trace(__VA_ARGS__)
	#define SAVANNAH_CORE_INFO(...)			::Savannah::Debugger::GetCoreLogger()->info(__VA_ARGS__)
	#define SAVANNAH_CORE_WARN(...)			::Savannah::Debugger::GetCoreDebugger()->warn(__VA_ARGS__)
	#define SAVANNAH_CORE_ERROR(...)		::Savannah::Debugger::GetCoreDebugger()->error(__VA_ARGS__)
	#define SAVANNAH_CORE_CRITICAL(...)		::Savannah::Debugger::GetCoreDebugger()->critical(__VA_ARGS__)
	
	// App log macros
	#define SAVANNAH_TRACE(...)		::Savannah::Debugger::GetClientDebugger()->trace(__VA_ARGS__)
	#define SAVANNAH_INFO(...)		::Savannah::Debugger::GetClientDebugger()->info(__VA_ARGS__)
	#define SAVANNAH_WARN(...)		::Savannah::Debugger::GetClientDebugger()->warn(__VA_ARGS__)
	#define SAVANNAH_ERROR(...)		::Savannah::Debugger::GetClientDebugger()->error(__VA_ARGS__)
	#define SAVANNAH_CRITICAL(...)	::Savannah::Debugger::GetClientDebugger()->critical(__VA_ARGS__)
#elif SAVANNAH_DEBUGGER
	// Savannah core log macros
	#define SAVANNAH_CORE_TRACE(...)		printf(__VA_ARGS__)
	#define SAVANNAH_CORE_INFO(...)			printf(__VA_ARGS__)
	#define SAVANNAH_CORE_WARN(...)			printf(__VA_ARGS__)
	#define SAVANNAH_CORE_ERROR(...)		printf(__VA_ARGS__)
	#define SAVANNAH_CORE_CRITICAL(...)		printf(__VA_ARGS__)
	
	// App log macros
	#define SAVANNAH_TRACE(...)		printf(__VA_ARGS__)
	#define SAVANNAH_INFO(...)		printf(__VA_ARGS__)
	#define SAVANNAH_WARN(...)		printf(__VA_ARGS__)
	#define SAVANNAH_ERROR(...)		printf(__VA_ARGS__)
	#define SAVANNAH_CRITICAL(...)	printf(__VA_ARGS__)
#else
	// Savannah core log macros
	#define SAVANNAH_CORE_TRACE(...)
	#define SAVANNAH_CORE_INFO(...)
	#define SAVANNAH_CORE_WARN(...)
	#define SAVANNAH_CORE_ERROR(...)
	#define SAVANNAH_CORE_CRITICAL(...)
	
	// App log macros
	#define SAVANNAH_TRACE(...)
	#define SAVANNAH_INFO(...)
	#define SAVANNAH_WARN(...)
	#define SAVANNAH_ERROR(...)
	#define SAVANNAH_CRITICAL(...)
#endif

#endif
