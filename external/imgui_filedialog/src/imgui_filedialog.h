/*
	Copyright 2020 Limeoats
	Original project: https://github.com/Limeoats/L2DFileDialog

	Changes by Vladimir Sigalkin
*/


#ifndef IMGUI_FILEDIALOG_H
#define IMGUI_FILEDIALOG_H

#include <filesystem>
#include <external/imgui/imgui.h>
#include <string>
#include <vector>

typedef int ImGuiFileDialogType;	// -> enum ImGuiFileDialogType_        // Enum: A file dialog type

enum ImGuiFileDialogType_
{
	ImGuiFileDialogType_OpenFile,
	ImGuiFileDialogType_SaveFile,
	ImGuiFileDialogType_COUNT
};

enum class ImGuiFileDialogSelectedType
{
	None,
	File,
	Directory
};

struct ImFileDialogInfo
{
	std::string title;
	ImGuiFileDialogType type;

	std::filesystem::path fileName;
	std::filesystem::path directoryPath;
	std::filesystem::path resultPath;
	std::vector<char*> fileExtensions = {};
	std::string fileExtension = "";
	int fileExtensionSelected = -1;
	ImGuiFileDialogSelectedType selection = ImGuiFileDialogSelectedType::None;

	bool refreshInfo;
	size_t currentIndex;
	std::vector<std::filesystem::directory_entry> currentFiles;
	std::vector<std::filesystem::directory_entry> currentDirectories;
//	std::function<bool(const std::string&)> fileFilterFunc = {};
	void RefreshPaths();
	bool FileFilter(const std::string& filename);
};

// ====================================================================
// for fileExtensions filter in ImFileDialogInfo::fileFilter()
#include <locale>  
#include <algorithm> 

template <typename T>
struct isEqualCaseInsensitive
{
	isEqualCaseInsensitive(const std::locale& loc) : m_locale(loc) {};
	bool operator()(T charA, T charB)
	{
		return std::toupper(charA, m_locale) == std::toupper(charB, m_locale);
	}
	
private:
	std::locale m_locale;
};

template <typename T>
int ci_findstr(const T& stringA, const T& stringB, const std::locale& locale = std::locale())
{
	typename T::const_iterator it = std::search(
		stringA.begin(), stringA.end(),
		stringB.begin(), stringB.end(),
		isEqualCaseInsensitive<typename T::value_type>(locale)
		);
	
	if (it != stringA.end())
	{
		return it - stringA.begin();
	}
	
	return -1; // nothing found
}
// ====================================================================

namespace ImGui
{
	IMGUI_API bool FileDialog(bool* open, ImFileDialogInfo* dialogInfo);
}

#endif
