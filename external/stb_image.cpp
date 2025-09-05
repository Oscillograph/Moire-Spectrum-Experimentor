#include <savannah/core.h>

#ifdef SAVANNAH_PLATFORM_WINDOWS
	// for stb_image.h : Unicode filenames on Windows
	#define STBI_WINDOWS_UTF8
#endif
#define STB_IMAGE_IMPLEMENTATION
#include <external/stb_image.h> // to work with images file formats

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <external/stb_image_write.h>
