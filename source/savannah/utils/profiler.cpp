#include <savannah/utils/profiler.h>

namespace Savannah {
	std::unordered_map<std::string, float> Profiler::m_ProfileResults;
	ProfilerSession* Profiler::m_CurrentSession = nullptr;
	std::string Profiler::m_Filepath = "";
	std::ofstream Profiler::m_OutputStream;
	int Profiler::m_ProfileCount = 0;
}
