/*
	Copyright 2020 Limeoats
	Original project: https://github.com/Limeoats/L2DFileDialog

	Changes by Vladimir Sigalkin
*/

/*
	Changes by Kirill Kamenskiy
	+ different languages support
	+ file dialog is now modal window instead of regular one
	+ mouse doubleclick on a file/directory processing
	+ Enter key press on selected file/directory processing
	+ file filtering reworked from user-defined lambda to predefined function
	+ file filtering is case-insensitive
	+ manual current path editing support
	+ various layout changes
	- date column is removed due to lack of knowledge about how to make it work
*/

#ifndef IMGUI_FILEDIALOG_CPP
#define IMGUI_FILEDIALOG_CPP

#include <chrono>
#include <string>
#include <filesystem>
#include <sstream>
#include <algorithm>
#include <cstring>

#include <external/imgui/imgui.h>
#include <external/imgui/misc/cpp/imgui_stdlib.h>
#include <external/imgui_filedialog/src/imgui_filedialog.h>
#include <external/imgui_filedialog/src/imgui_filedialog_lang.h>

using namespace std::chrono_literals;

#include <logger.h> // Savannah App debugging

typedef int ImGuiFileDialogSortOrder;

enum ImGuiFileDialogSortOrder_
{
	ImGuiFileDialogSortOrder_Up,
	ImGuiFileDialogSortOrder_Down,
	ImGuiFileDialogSortOrder_None
};

void ImFileDialogInfo::RefreshPaths()
{
	refreshInfo = false;
	currentIndex = 0;
	currentFiles.clear();
	currentDirectories.clear();

	for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(directoryPath))
	{
		if (entry.is_directory())
		{
			currentDirectories.push_back(entry);
		}
		else
		{
			if (fileExtensionSelected != -1) 
			{
				if (!FileFilter(entry.path().filename().string())) 
				{
					continue;
				}
			}
			currentFiles.push_back(entry);
		}
	}
}

bool ImFileDialogInfo::FileFilter(const std::string& filename)
{
	bool result = false;
	std::size_t found;
	if (fileExtensionSelected != -1)
	{
		found = ci_findstr(filename, std::string(fileExtensions[fileExtensionSelected]));
		if (found != std::string::npos)
		{
			result = true;
//			break;
		}
	} else {
		result = true;
	}
	return result;
}

bool ImGui::FileDialog(bool* open, ImFileDialogInfo* dialogInfo)
{
	if (!*open) return false;

	static float initialSpacingColumn0 = 430.0f;
	static float initialSpacingColumn1 = 80.0f;
	static float initialSpacingColumn2 = 90.0f;
	static ImGuiFileDialogSortOrder fileNameSortOrder = ImGuiFileDialogSortOrder_None;
	static ImGuiFileDialogSortOrder sizeSortOrder = ImGuiFileDialogSortOrder_None;
	static ImGuiFileDialogSortOrder dateSortOrder = ImGuiFileDialogSortOrder_None;
	static ImGuiFileDialogSortOrder typeSortOrder = ImGuiFileDialogSortOrder_None;
	
	static std::string currentDirectoryPathString = dialogInfo->directoryPath.string();

	assert(dialogInfo != nullptr);

	bool complete = false;
	bool fileDoubleClick = false;
	bool directoryDoubleClick = false;

	if (!ImGui::fileDialogLanguagesInitialized)
	{
		ImGui::FileDialogLanguageInitialize();
	}
	
	ImGui::PushID(dialogInfo);
	ImGui::SetNextWindowSize(ImVec2(740.0f, 420.0f), ImGuiCond_FirstUseEver);

	ImGui::OpenPopup(dialogInfo->title.c_str());
	if (ImGui::BeginPopupModal(dialogInfo->title.c_str(), open, ImGuiWindowFlags_NoResize))
	{
		if (dialogInfo->currentFiles.empty() && dialogInfo->currentDirectories.empty() || dialogInfo->refreshInfo)
			dialogInfo->RefreshPaths();
		
		// handle pressed Enter button event
		if (ImGui::IsKeyDown(ImGuiKey::ImGuiKey_Enter))
		{
			if (dialogInfo->selection == ImGuiFileDialogSelectedType::File)
			{
				fileDoubleClick = true;
			}
			if (dialogInfo->selection == ImGuiFileDialogSelectedType::Directory)
			{
				directoryDoubleClick = true;
			}
		}

		// Draw path
//		ImGui::Text("%s: %s", ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Path").c_str(), dialogInfo->directoryPath.string().c_str());
		ImGui::Text("%s:", ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Path").c_str());
		ImGui::SameLine();
		ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
		if (ImGui::InputText("###InputTextCurrentPath", &(currentDirectoryPathString), 4096) || ImGui::IsKeyDown(ImGuiKey::ImGuiKey_Enter))
		{
			std::filesystem::path path(currentDirectoryPathString);
			if (std::filesystem::exists(path))
			{
				dialogInfo->directoryPath = path;
				dialogInfo->RefreshPaths();
			}
		}

		ImGui::BeginChild("##browser", ImVec2(ImGui::GetWindowContentRegionWidth(), 290), true, ImGuiWindowFlags_HorizontalScrollbar);
		ImGui::Columns(3);

		// Columns size
		if (initialSpacingColumn0 > 0)
		{
			ImGui::SetColumnWidth(0, initialSpacingColumn0);
			initialSpacingColumn0 = 0.0f;
		}
		if (initialSpacingColumn1 > 0)
		{
			ImGui::SetColumnWidth(1, initialSpacingColumn1);
			initialSpacingColumn1 = 0.0f;
		}
//		if (initialSpacingColumn2 > 0)
//		{
//			ImGui::SetColumnWidth(2, initialSpacingColumn2);
//			initialSpacingColumn2 = 0.0f;
//		}

		// File Columns
		if (ImGui::Selectable(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Name").c_str()))
		{
			sizeSortOrder = ImGuiFileDialogSortOrder_None;
			dateSortOrder = ImGuiFileDialogSortOrder_None;
			typeSortOrder = ImGuiFileDialogSortOrder_None;
			fileNameSortOrder = fileNameSortOrder == ImGuiFileDialogSortOrder_Down ? ImGuiFileDialogSortOrder_Up : ImGuiFileDialogSortOrder_Down;
		}
		ImGui::NextColumn();
		if (ImGui::Selectable(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Size").c_str()))
		{
			fileNameSortOrder = ImGuiFileDialogSortOrder_None;
			dateSortOrder = ImGuiFileDialogSortOrder_None;
			typeSortOrder = ImGuiFileDialogSortOrder_None;
			sizeSortOrder = sizeSortOrder == ImGuiFileDialogSortOrder_Down ? ImGuiFileDialogSortOrder_Up : ImGuiFileDialogSortOrder_Down;
		}
		ImGui::NextColumn();
		if (ImGui::Selectable(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Type").c_str()))
		{
			fileNameSortOrder = ImGuiFileDialogSortOrder_None;
			dateSortOrder = ImGuiFileDialogSortOrder_None;
			sizeSortOrder = ImGuiFileDialogSortOrder_None;
			typeSortOrder = typeSortOrder == ImGuiFileDialogSortOrder_Down ? ImGuiFileDialogSortOrder_Up : ImGuiFileDialogSortOrder_Down;
		}
		ImGui::NextColumn();
//		if (ImGui::Selectable(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Date").c_str()))
//		{
//			fileNameSortOrder = ImGuiFileDialogSortOrder_None;
//			sizeSortOrder = ImGuiFileDialogSortOrder_None;
//			typeSortOrder = ImGuiFileDialogSortOrder_None;
//			dateSortOrder = dateSortOrder == ImGuiFileDialogSortOrder_Down ? ImGuiFileDialogSortOrder_Up : ImGuiFileDialogSortOrder_Down;
//		}
//		ImGui::NextColumn();

		// File Separator
		ImGui::Separator();

		// Sort directories
		auto* directories = &dialogInfo->currentDirectories;

		if (fileNameSortOrder != ImGuiFileDialogSortOrder_None || sizeSortOrder != ImGuiFileDialogSortOrder_None || typeSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(directories->begin(), directories->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (fileNameSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.path().filename() > b.path().filename();
				}

				return a.path().filename() < b.path().filename();
			});
		}
		else if (dateSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(directories->begin(), directories->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (dateSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.last_write_time() > b.last_write_time();
				}
				return a.last_write_time() < b.last_write_time();
			});
		}

		// Sort files
		auto* files = &dialogInfo->currentFiles;

		if (fileNameSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(files->begin(), files->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (fileNameSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.path().filename() > b.path().filename();
				}

				return a.path().filename() < b.path().filename();
			});
		}
		else if (sizeSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(files->begin(), files->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (sizeSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.file_size() > b.file_size();
				}
				return a.file_size() < b.file_size();
			});
		}
		else if (typeSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(files->begin(), files->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (typeSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.path().extension() > b.path().extension();
				}

				return a.path().extension() < b.path().extension();
			});
		}
		else if (dateSortOrder != ImGuiFileDialogSortOrder_None)
		{
			std::sort(files->begin(), files->end(), [](const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b)
			{
				if (dateSortOrder == ImGuiFileDialogSortOrder_Down)
				{
					return a.last_write_time() > b.last_write_time();
				}
				return a.last_write_time() < b.last_write_time();
			});
		}

		size_t index = 0;

		// Draw parent
		if (dialogInfo->directoryPath.has_parent_path())
		{
			if (ImGui::Selectable("..", dialogInfo->currentIndex == index, ImGuiSelectableFlags_AllowDoubleClick, ImVec2(ImGui::GetWindowContentRegionWidth(), 0)))
			{
				dialogInfo->currentIndex = index;
				dialogInfo->selection = ImGuiFileDialogSelectedType::Directory;

				if (ImGui::IsMouseDoubleClicked(0) || directoryDoubleClick)
				{
					dialogInfo->directoryPath = dialogInfo->directoryPath.parent_path();
					currentDirectoryPathString = dialogInfo->directoryPath.string();
					dialogInfo->refreshInfo = true;
					dialogInfo->selection = ImGuiFileDialogSelectedType::None;
				}
			}
			
			ImGui::NextColumn();
			ImGui::TextUnformatted("-");
			ImGui::NextColumn();
			ImGui::TextUnformatted(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "<parent>").c_str());
			ImGui::NextColumn();
//			ImGui::TextUnformatted("-");
//			ImGui::NextColumn();

			index++;
		}

		// Draw directories
		for (size_t i = 0; i < directories->size(); ++i)
		{
			auto directoryEntry = dialogInfo->currentDirectories[i];
			auto directoryPath = directoryEntry.path();
			auto directoryName = directoryPath.filename();

			if (ImGui::Selectable(directoryName.string().c_str(), dialogInfo->currentIndex == index, ImGuiSelectableFlags_AllowDoubleClick, ImVec2(ImGui::GetWindowContentRegionWidth(), 0)))
			{
				dialogInfo->currentIndex = index;
				dialogInfo->selection = ImGuiFileDialogSelectedType::Directory;

				if (ImGui::IsMouseDoubleClicked(0) || directoryDoubleClick)
				{
					dialogInfo->directoryPath = directoryPath;
					currentDirectoryPathString = dialogInfo->directoryPath.string();
					dialogInfo->refreshInfo = true;
					dialogInfo->selection = ImGuiFileDialogSelectedType::None;
				}
			}

			ImGui::NextColumn();
			ImGui::TextUnformatted("-");
			ImGui::NextColumn();
			ImGui::TextUnformatted(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "<directory>").c_str());
			ImGui::NextColumn();

//			std::filesystem::file_time_type lastWriteTime = std::filesystem::last_write_time(directoryPath);// directoryEntry.last_write_time();
//			auto st = std::chrono::time_point_cast<std::chrono::system_clock::duration>(lastWriteTime - decltype(lastWriteTime)::clock::now() + std::chrono::system_clock::now());
//			std::time_t tt = std::chrono::system_clock::to_time_t(st);
//			std::tm* mt;
//
//			#ifdef _WIN32
//			localtime_s(mt, &tt);
//			#else
//			mt = localtime(&tt);
//			#endif
//
//			std::stringstream ss;
//			ss << std::put_time(mt, "%F %R");
//
//			ImGui::TextUnformatted(ss.str().c_str());
//			ImGui::NextColumn();

			index++;
		}

		// Draw files
		for (size_t i = 0; i < files->size(); ++i)
		{
			auto fileEntry = dialogInfo->currentFiles[i];
			auto filePath = fileEntry.path();
			auto fileName = filePath.filename();

			if (ImGui::Selectable(fileName.string().c_str(), dialogInfo->currentIndex == index, ImGuiSelectableFlags_AllowDoubleClick, ImVec2(ImGui::GetWindowContentRegionWidth(), 0)))
			{
				dialogInfo->currentIndex = index;
				dialogInfo->selection = ImGuiFileDialogSelectedType::File;
				dialogInfo->fileName = fileName;
				
				if (ImGui::IsMouseDoubleClicked(0))
				{
					fileDoubleClick = true;
				}
			}

			ImGui::NextColumn();
			ImGui::TextUnformatted(std::to_string(fileEntry.file_size()).c_str());
			ImGui::NextColumn();
			ImGui::TextUnformatted(filePath.extension().string().c_str());
			ImGui::NextColumn();

//			auto lastWriteTime = fileEntry.last_write_time();
//			auto st = std::chrono::time_point_cast<std::chrono::system_clock::duration>(lastWriteTime - decltype(lastWriteTime)::clock::now() + std::chrono::system_clock::now());
//			std::time_t tt = std::chrono::system_clock::to_time_t(st);
//			std::tm* mt;
//			
//			#ifdef _WIN32
//			localtime_s(mt, &tt);
//			#else
//			mt = localtime(&tt);
//			#endif
//
//			std::stringstream ss;
//			ss << std::put_time(mt, "%F %R");
//			ImGui::TextUnformatted(ss.str().c_str());
//			ImGui::NextColumn();

			index++;
		}
		
		ImGui::EndChild();

		// Draw filename
		static const size_t fileNameBufferSize = 200;
		static char fileNameBuffer[fileNameBufferSize];

		std::string fileNameStr = dialogInfo->fileName.string();
		size_t fileNameSize = fileNameStr.size();

		if (fileNameSize >= fileNameBufferSize)	fileNameSize = fileNameBufferSize - 1;
		std::memcpy(fileNameBuffer, fileNameStr.c_str(), fileNameSize);
		fileNameBuffer[fileNameSize] = 0;

		ImGui::PushItemWidth(150);
		ImGui::Text(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "File Name").c_str());
		ImGui::SameLine();
		
		ImGui::PushItemWidth(ImGui::GetWindowContentRegionWidth() - 250);
		if (ImGui::InputText("###InputTextFileName", fileNameBuffer, fileNameBufferSize))
		{
			dialogInfo->fileName = std::string(fileNameBuffer);
			dialogInfo->currentIndex = 0;
		}
		ImGui::SameLine();

		std::string fileTypePreview;
		if (dialogInfo->fileExtensionSelected != -1)
		{
			fileTypePreview = dialogInfo->fileExtensions[dialogInfo->fileExtensionSelected];
		} else {
			fileTypePreview = ".*";
		}
		
		ImGui::PushItemWidth(150);
		if (ImGui::BeginCombo("###FileType", fileTypePreview.c_str()))
		{
			bool isSelected = false;
			
			for (int i = -1; i < (int)(dialogInfo->fileExtensions.size()); ++i)
			{
				isSelected = (dialogInfo->fileExtensionSelected == i);
				if (i == -1)
				{
					if (ImGui::Selectable(".*", isSelected))
					{
						dialogInfo->fileExtensionSelected = -1;
						dialogInfo->RefreshPaths();
					}
				} else {
					if (ImGui::Selectable(dialogInfo->fileExtensions[i], isSelected))
					{
						dialogInfo->fileExtensionSelected = i;
						dialogInfo->RefreshPaths();
					}
				}
				
				if (isSelected)
				{
					ImGui::SetItemDefaultFocus();
				}
			}
			
			ImGui::EndCombo();
		}
		

//		ImGui::SetNextItemWidth();
//		ImGui::Separator();
		ImGui::NewLine();
		ImGui::SameLine(ImGui::GetContentRegionAvail().x - 150);
		ImGui::SetNextItemWidth(75);
		if (ImGui::Button(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Cancel").c_str(), {70, 25}))
		{
			fileNameSortOrder = ImGuiFileDialogSortOrder_None;
			sizeSortOrder = ImGuiFileDialogSortOrder_None;
			typeSortOrder = ImGuiFileDialogSortOrder_None;
			dateSortOrder = ImGuiFileDialogSortOrder_None;

			dialogInfo->refreshInfo = false;
			dialogInfo->currentIndex = 0;
			dialogInfo->currentFiles.clear();
			dialogInfo->currentDirectories.clear();

			*open = false;
			ImGui::CloseCurrentPopup();
		}

		ImGui::SameLine(ImGui::GetContentRegionAvail().x - 70);
		ImGui::SetNextItemWidth(75);
		if (dialogInfo->type == ImGuiFileDialogType_OpenFile)
		{
			if (ImGui::Button(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Open").c_str(), {70, 25}) || fileDoubleClick)
			{
				dialogInfo->resultPath = dialogInfo->directoryPath / dialogInfo->fileName;

				if (std::filesystem::exists(dialogInfo->resultPath))
				{
					fileNameSortOrder = ImGuiFileDialogSortOrder_None;
					sizeSortOrder = ImGuiFileDialogSortOrder_None;
					typeSortOrder = ImGuiFileDialogSortOrder_None;
					dateSortOrder = ImGuiFileDialogSortOrder_None;

					dialogInfo->refreshInfo = false;
					dialogInfo->currentIndex = 0;
					dialogInfo->currentFiles.clear();
					dialogInfo->currentDirectories.clear();
					dialogInfo->fileExtension = dialogInfo->resultPath.extension().string();

					complete = true;
					*open = false;
					ImGui::CloseCurrentPopup();
				}
			}
		}
		else if (dialogInfo->type == ImGuiFileDialogType_SaveFile)
		{
			if (ImGui::Button(ImGui::FileDialogLanguageGetText(ImGui::fileDialogLanguageSelected, "Save").c_str()))
			{
				dialogInfo->resultPath = dialogInfo->directoryPath / dialogInfo->fileName;

				if (!std::filesystem::exists(dialogInfo->resultPath))
				{
					fileNameSortOrder = ImGuiFileDialogSortOrder_None;
					sizeSortOrder = ImGuiFileDialogSortOrder_None;
					typeSortOrder = ImGuiFileDialogSortOrder_None;
					dateSortOrder = ImGuiFileDialogSortOrder_None;

					dialogInfo->refreshInfo = false;
					dialogInfo->currentIndex = 0;
					dialogInfo->currentFiles.clear();
					dialogInfo->currentDirectories.clear();
					dialogInfo->fileExtension = dialogInfo->resultPath.extension().string();

					complete = true;
					*open = false;
					ImGui::CloseCurrentPopup();
				}
			}
		}
		
		ImGui::EndPopup();
	}

	ImGui::PopID();

	return complete;
}

#endif
