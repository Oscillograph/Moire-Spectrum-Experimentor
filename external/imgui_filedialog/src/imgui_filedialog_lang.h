#ifndef IMGUI_FILEDIALOG_LANG_H
#define IMGUI_FILEDIALOG_LANG_H

#include <string>
#include <map>

namespace ImGui
{
	enum class FileDialogLanguage
	{
		None				,
		English				,
		Russian				,
	};
	
	typedef std::map<std::string, std::string> LanguageMap;
	bool fileDialogLanguagesInitialized = false;
	FileDialogLanguage fileDialogLanguageSelected = FileDialogLanguage::None;
	std::map<FileDialogLanguage, LanguageMap> FileDialogLanguageText = {};
	
	void FileDialogLanguageInitialize(FileDialogLanguage lang = FileDialogLanguage::Russian)
	{
		// Default
		LanguageMap noneLanguageMap =
		{
			{"no text", 		"%%%"},
		};
		FileDialogLanguageText[FileDialogLanguage::None] = noneLanguageMap;
		
		// Russian
		LanguageMap russianLanguageMap = 
		{
			{"no text", 		"нет текста"},
			{"Path", 			"Путь"},
			{"Name", 			"Имя"},
			{"Size", 			"Размер"},
			{"Type", 			"Тип"},
			{"<parent>",		"<вверх>"},
			{"<directory>",		"<папка>"},
			{"Date", 			"Дата последнего изменения"},
			{"File Name", 		"Имя файла"},
			{"Cancel", 			"Отмена"},
			{"Open", 			"Открыть"},
			{"Save", 			"Сохранить"},
		};
		FileDialogLanguageText[FileDialogLanguage::Russian] = russianLanguageMap;
		
		// English
		LanguageMap englishLanguageMap =
		{
			{"no text", 		"no text"},
			{"Path", 			"Path"},
			{"Name", 			"Name"},
			{"Size", 			"Size"},
			{"Type", 			"Type"},
			{"<parent>",		"<parent>"},
			{"<directory>",		"<directory>"},
			{"Date", 			"Date"},
			{"File Name", 		"File Name"},
			{"Cancel", 			"Cancel"},
			{"Open", 			"Open"},
			{"Save", 			"Save"},
		};
		
		FileDialogLanguageText[FileDialogLanguage::English] = englishLanguageMap;
		
		fileDialogLanguageSelected = lang;
		fileDialogLanguagesInitialized = true;
	}
	
	std::string& FileDialogLanguageGetText(FileDialogLanguage lang, const std::string& name)
	{
		if (FileDialogLanguageText.find(lang) != FileDialogLanguageText.end())
		{
			if (FileDialogLanguageText[lang].find(name) != FileDialogLanguageText[lang].end())
			{
				return FileDialogLanguageText[lang][name];
			} else {
				return FileDialogLanguageText[FileDialogLanguage::None]["no-text"];
			}
		} else {
			return FileDialogLanguageText[FileDialogLanguage::None]["no-text"];
		}
	}
}

#endif
