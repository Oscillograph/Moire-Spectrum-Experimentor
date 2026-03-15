#ifndef SAVANNAH_LOGGER_H
#define SAVANNAH_LOGGER_H

#include <savannah/proto-core.h>
#include <savannah/utils/time.h>
#include <external/utf8.h>

// comment this line if you want something simple like std::cout
// #define SAVANNAH_LOGGER_CLASS

namespace Savannah
{
	typedef std::basic_stringstream<char32_t> u32stringstream;
	
	class Logger
	{
	public:
		static void Flush(int mode = 0, int color = -1);
		
		template<typename... Args>
		inline static void Text(Args... args)
		{
			m_u8stream.str(std::string()); // clear utf8 buffer
			m_u32string = U""; // clear utf32 string
			
			// add args to utf8 string stream
			m_u8stream << "[" << Savannah::GetDateAndTime() << "] ";
			(u8add(std::forward<Args>(args)), ...);
			m_u8stream << "\n";
			m_u8string = m_u8stream.str();
			
			// convert it to utf32
			utf8::utf8to32(m_u8string.begin(), m_u8string.end(), std::back_inserter(m_u32string));
			
			m_u32stream << m_u32string << "\n";
		}
		
		template<typename... Args>
		inline static void Message(Args... args)
		{
			m_u8stream.str(std::string()); // clear utf8 buffer
			m_u32string = U""; // clear utf32 string
			
			// add args to utf8 string stream
			// m_u8stream << "[" << Utils::GetDateAndTime() << "] ";
			(u8add(std::forward<Args>(args)), ...);
			m_u8stream << "\n";
			m_u8string = m_u8stream.str();
			
			// convert it to utf32
			utf8::utf8to32(m_u8string.begin(), m_u8string.end(), std::back_inserter(m_u32string));
			
			m_u32stream << m_u32string << "\n";
		}
		
	private:
		inline static void Reset()
		{
			m_u8string = "";
			m_u32string = U"";
			m_u8stream.str(std::string());
			m_u32stream.str(std::u32string());
		}
		
		template <typename T>
		inline static void u8add(T t) { m_u8stream << t; }
		
		template <typename T>
		inline static void u32add(T t) { m_u32stream << t; }
		
		static std::stringstream m_u8stream;
		static std::string m_u8string;
		static std::u32string m_u32string;
		static u32stringstream m_u32stream;
	};
}

#ifndef SAVANNAH_LOGGER_CLASS
	#define SAVANNAH_CONSOLE_LOG(...) 				std::cout << "App Log: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CONSOLE_DEBUG(...)				std::cout << "App Debug: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CONSOLE_ERROR(...)				std::cout << "App Error: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CONSOLE_ASSERT(x, ...)			{ if (!(x)) { std::cout << "App Assert: " << __VA_ARGS__ (__VA_ARGS__); std::exit(); } }
	
	#define SAVANNAH_CORE_TRACE(...)				std::cout << "Core Trace: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CORE_INFO(...)					std::cout << "Core Info: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CORE_WARN(...)					std::cout << "Core Warning: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CORE_ERROR(...)				std::cout << "Core Error: " << __VA_ARGS__ << "\n";
	#define SAVANNAH_CORE_CRITICAL(...)				std::cout << "Core Critical: " << __VA_ARGS__ << "\n";
#else
	// Savannah misc console macros
	#define SAVANNAH_CONSOLE_GREY(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, -1);
	#define SAVANNAH_CONSOLE_YELLOW(...)			::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 0);
	#define SAVANNAH_CONSOLE_WHITE(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 1);
	#define SAVANNAH_CONSOLE_DARKRED(...)			::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 2);
	#define SAVANNAH_CONSOLE_RED(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 3);
	#define SAVANNAH_CONSOLE_TEAL(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 4);
	#define SAVANNAH_CONSOLE_CYAN(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 5);
	#define SAVANNAH_CONSOLE_BLUE(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 6);
	#define SAVANNAH_CONSOLE_MAGENTA(...)			::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 7);
	#define SAVANNAH_CONSOLE_GREEN(...)				::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 8);
	#define SAVANNAH_CONSOLE_CAPTION_BLACK(...)		::Savannah::Logger::Message(__VA_ARGS__); ::Savannah::Logger::Flush(0, 9);
	
	#define SAVANNAH_CONSOLE_LOG(...)				::Savannah::Logger::Text("App: ", __VA_ARGS__); ::Savannah::Logger::Flush(0, 1);
	#define SAVANNAH_CONSOLE_ERROR(...)				::Savannah::Logger::Text("App: ", __VA_ARGS__); ::Savannah::Logger::Flush(0, 3);
	#define SAVANNAH_CONSOLE_ASSERT(x, ...)			{ if (!(x)) { SAVANNAH_CORE_ERROR(__VA_ARGS__); std::exit(-4); } }
	#define SAVANNAH_CONSOLE_DEBUG(...)				::Savannah::Logger::Text("App: ", __VA_ARGS__); ::Savannah::Logger::Flush(1, 1);

	// Savannah core log macros
	#define SAVANNAH_CORE_TRACE(...)		SAVANNAH_CONSOLE_GREY("Core Trace: ", __VA_ARGS__)
	#define SAVANNAH_CORE_INFO(...)			SAVANNAH_CONSOLE_GREY("Core Info: ", __VA_ARGS__)
	#define SAVANNAH_CORE_WARN(...)			SAVANNAH_CONSOLE_YELLOW("Core Warning: ", __VA_ARGS__)
	#define SAVANNAH_CORE_ERROR(...)		SAVANNAH_CONSOLE_RED("Core Error: ", __VA_ARGS__)
	#define SAVANNAH_CORE_CRITICAL(...)		SAVANNAH_CONSOLE_RED("Core Critical Error: ", __VA_ARGS__)
	
	// App log macros
	#define SAVANNAH_TRACE(...)		SAVANNAH_CONSOLE_GREY("App Trace: ", __VA_ARGS__)
	#define SAVANNAH_INFO(...)		SAVANNAH_CONSOLE_GREY("App Info: ", __VA_ARGS__)
	#define SAVANNAH_WARN(...)		SAVANNAH_CONSOLE_YELLOW("App Warning: ", __VA_ARGS__)
	#define SAVANNAH_ERROR(...)		SAVANNAH_CONSOLE_RED("App Error: ", __VA_ARGS__)
	#define SAVANNAH_CRITICAL(...)	SAVANNAH_CONSOLE_RED("App Critical Error: ", __VA_ARGS__)
#endif

#endif

