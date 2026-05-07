#include <savannah/shader.h>
#include <savannah/platforms/opengl/opengl_shader.h>

#include <savannah/utils/file-io.h>

namespace Savannah
{

	// ----------------------------------------------------------------
	//								Shader
	// ----------------------------------------------------------------
	Shader* Shader::Create(const std::string& name, const std::string& vertexSrc, const std::string& fragmentSrc){
		return new OpenGLShader(name, vertexSrc, fragmentSrc);
	}

	Shader* Shader::Create(const std::string& name, const std::string& filepath){
		return new OpenGLShader(name, filepath);
	}


	// ----------------------------------------------------------------
	//							Shader Library
	// ----------------------------------------------------------------
	void ShaderLibrary::Add(const std::string& name, Shader* shader){
		SAVANNAH_CORE_ASSERT(m_Shaders.find(name) == m_Shaders.end(), "Shader already exists in the library");
		m_Shaders[name] = shader;
	}

	void ShaderLibrary::Add(Shader* shader){
		auto& name = shader->GetName();
		SAVANNAH_CORE_ASSERT(m_Shaders.find(name) == m_Shaders.end(), "Shader already exists in the library");
		m_Shaders[name] = shader;
	}

	Shader* ShaderLibrary::Load(const std::string& filepath){
		auto name = FileIO::GetName(filepath);
		auto shader = Shader::Create(name, filepath);
		Add(shader);
		return m_Shaders[name];
	}

	Shader* ShaderLibrary::Load(const std::string& name, const std::string& filepath){
		auto shader = Shader::Create(name, filepath);
		Add(name, shader);
		return m_Shaders[name];
	}

	Shader* ShaderLibrary::Get(const std::string& name){
		if (m_Shaders.find(name) != m_Shaders.end()){
			return m_Shaders[name];
		} else {
			return nullptr;
		}
		SAVANNAH_CORE_ASSERT(m_Shaders.find(name) != m_Shaders.end(), "Shader doesn't exist in the library");
	}
}
