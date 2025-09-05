#ifndef SAVANNAH_PROTOCORE_H
#define SAVANNAH_PROTOCORE_H

// -------------------------------------------------------------------------------
//     The point of this "Proto" is to fix the multiple definition error that
// occurs through log.cpp file, which includes Core.h while Core.h includes log.h
// -------------------------------------------------------------------------------

#ifdef SAVANNAH_PLATFORM_WINDOWS
	#define SAVANNAH_API
#else
	#error Savannah currently only supports Windows.
#endif

// Assertions for debugging
#ifdef SAVANNAH_ENABLE_ASSERTS
	#define SAVANNAH_ASSERT(x,...) 		{if(!(x)) SAVANNAH_ERROR("Assertion failed: {0}", __VA_ARGS__); __debugbreak(); }
	#define SAVANNAH_CORE_ASSERT(x,...) {if(!(x)) SAVANNAH_CORE_ERROR("Assertion failed: {0}", __VA_ARGS__); __debugbreak(); }
#else
	#define SAVANNAH_ASSERT(x,...)
	#define SAVANNAH_CORE_ASSERT(x,...)
#endif

// useful macros
#define BIT(x) (1 << x)

// FPS							ms between frames 
#define SAVANNAH_FPSNOLIMIT		0.0f
#define SAVANNAH_FPS120			0.0083f
#define SAVANNAH_FPS90			0.0111f
#define SAVANNAH_FPS60			0.0167f
#define SAVANNAH_FPS30			0.0333f
#define SAVANNAH_FPS15			0.0667f
#define SAVANNAH_FPS10			0.1f
#define SAVANNAH_FPS5			0.2f
#define SAVANNAH_FPS4			0.25f
#define SAVANNAH_FPS3			0.3333f
#define SAVANNAH_FPS2			0.5f
#define SAVANNAH_FPS1			1.0f

// defines before anything Vendor-ish is included
// #define GLM_FORCE_CTOR_INIT // to force GLM initialize vectors and matrices 

// basic includes
#include <savannah/common.h>

namespace Savannah {
	// Still no idea why do I need this typedef-kind-of stuff.
	// The Cherno adds it, but I don't get the point.
	template <typename T>
	using Scope = std::unique_ptr<T>;
	
	template <typename T>
	using Ref = std::shared_ptr<T>;
	// Thing is, shared pointers make use of a reference counter which is incremented through
	// operator "=" and std::make_shared<T>(value) function.
	// Unique pointers are created with constructor std::unique_ptr<T> m(new T(value)).
	// They have no reference count, and they get destroyed when they are out of scope.
	// So, usefulness of these aliases is doubtful at best. 
}

// basic types
namespace Savannah
{
	struct iVec2
	{
		int x = 0;
		int y = 0;
	};
	
	struct uVec2
	{
		unsigned int x = 0;
		unsigned int y = 0;
	};
	
	struct fVec2
	{
		float x = 0.0f;
		float y = 0.0f;
	};
	
	struct dVec2
	{
		double x = 0.0;
		double y = 0.0;
	};
}

#endif
