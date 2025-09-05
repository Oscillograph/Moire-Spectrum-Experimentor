#ifndef SAVANNAH_GAMEDATA_H
#define SAVANNAH_GAMEDATA_H

// =================================================================================================
// EngineData is a global static class to provide API for BSE and Client App systems to 
// communicate with each other despite architecture restrictions.
// It is a sort of a backdoor for circumstances where "You can't, but you you can if you want".
//
// To be able to store any object this class should be defined after everything else.
// Sadly, it can't store Application stuff, but that is to be defined in a client descendant anyway.
// =================================================================================================

#include <savannah/core.h>

namespace Savannah {
	struct SAVANNAH_API EngineData {
		// these should be copied from Application
		// static Window* m_Window;
		// static LayerStack* m_LayerStack;
		// static Layer* m_ImGuiLayer;
		// static FrameBuffer* m_FrameBuffer;
		// static bool m_ImGuiLayerEnabled;
		static float FPS_deltaTime;
		// static glm::vec4 RendererClearColor;
		// static glm::vec4 CustomColor;
		// static OrthographicCameraController* m_CameraController;
	};
}

#endif
