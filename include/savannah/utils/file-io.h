#ifndef SAVANNAH_FILEIO_H
#define SAVANNAH_FILEIO_H

#include <savannah/core.h>

namespace Savannah {
	class SAVANNAH_API FileIO {
	public:
		//static void ConstructFileSystem();
		
		static std::string GetRawText(std::string filename);
		
		static bool WriteRawText(const std::string& filename, const char* text);
		
		static std::string GetName(std::string filepath);
	};
}

#endif
