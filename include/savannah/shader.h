#ifndef SAVANNAH_SHADER
#define SAVANNAH_SHADER

#include <savannah/proto-core.h>

#include <external/glm/glm.hpp> // for vec and mat template types

namespace Savannah
{
	class Shader
	{
	public:
		virtual ~Shader() {};

		virtual void Bind() const = 0;
		virtual void Unbind() const = 0;

		inline void SetRendererID(uint32_t rendererId) { m_RendererId = rendererId; }
		inline uint32_t GetRendererID() const { return m_RendererId; }

		// virtual void UploadUniformBuffer(std::string name, T buffer) = 0;
		// virtual void UploadUniformBuffer() = 0;

		virtual void UploadUniformInt(std::string name, const int& value) const = 0;
		virtual void UploadUniformIntArray(std::string name, int* values, int count) const = 0;
		virtual void UploadUniformFloat(std::string name, const float& value) const = 0;
		virtual void UploadUniformFloatArray(std::string name, float* values, int count) const = 0;

		virtual void UploadUniformVec2(std::string name, const glm::vec2& vector) const = 0;
		virtual void UploadUniformVec3(std::string name, const glm::vec3& vector) const = 0;
		virtual void UploadUniformVec4(std::string name, const glm::vec4& vector) const = 0;

		virtual void UploadUniformMat3(std::string name, const glm::mat3& matrix) const = 0;
		virtual void UploadUniformMat4(std::string name, const glm::mat4& matrix) const = 0;

		inline const std::string& GetName() const { return m_Name; }
		inline void SetName(std::string name) { m_Name = name; }

		static Shader* Create(const std::string& name, const std::string& vertexSrc, const std::string& fragmentSrc);
		static Shader* Create(const std::string& name, const std::string& filepath); // for glsl-files with two shader programs in'em

	protected:
		uint32_t m_RendererId;
		std::string m_Name;
	};


	class ShaderLibrary
	{
	public:
		void Add(const std::string& name, Shader* shader);
		void Add(Shader* shader);
		Shader* Load(const std::string& filepath); // assets/shaders/(. . .).glsl
		Shader* Load(const std::string& name, const std::string& filepath);
		Shader* Get(const std::string& name); // geta shader from a library by name

	private:
		std::unordered_map<std::string, Shader*> m_Shaders;
	};
}

#endif
