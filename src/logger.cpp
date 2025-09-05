#include "../include/logger.h"
#include <external/color-console/color.hpp>

std::string Logger::m_u8string = "";
std::u32string Logger::m_u32string = U"";
std::stringstream Logger::m_u8stream;
u32stringstream Logger::m_u32stream;

void Logger::Flush(int mode, int color){
	// flush into command line interface
	if (mode == 0) 
	{
		switch (color)
		{
			case 0: // engine log
			{
				std::cout << dye::light_yellow(m_u8string);
				break;
			} 
			
			case 1: // app log
			{
				std::cout << dye::white(m_u8string);
				break;
			}
			
			case 2: // dark red
			{
				std::cout << dye::red(m_u8string);
				break;
			}
			
			case 3: // light red
			{
				std::cout << dye::light_red(m_u8string);
				break;
			}
			
			case 4: // teal
			{
				std::cout << dye::aqua(m_u8string);
				break;
			}
			
			case 5: // cyan
			{
				std::cout << dye::light_aqua(m_u8string);
				break;
			}
			
			case 6: // blue
			{
				std::cout << dye::light_blue(m_u8string);
				break;
			}
			
			case 7: // magenta
			{
				std::cout << dye::light_purple(m_u8string);
				break;
			}
			
			case 8: // green
			{
				std::cout << dye::light_green(m_u8string);
				break;
			}
			
			case 9: // caption: black on white
			{
				std::cout << dye::black_on_white(m_u8string);
				break;
			}
			
			default: // no idea who's log
			{
				std::cout << dye::grey(m_u8string);
			}
		}
	}
	
	// flush onto screen
	if (mode == 1) {}
	
	Reset();
}
