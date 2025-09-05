#ifndef SAVANNAH_PLATFORMUTILS_H
#define SAVANNAH_PLATFORMUTILS_H

#include <savannah/proto-core.h>

namespace BSE {
	class SAVANNAH_API FileDialogs {
	public:
		static std::string OpenFile(const char* filter);
		static std::string SaveFile(const char* filter);
	};
}

#endif
