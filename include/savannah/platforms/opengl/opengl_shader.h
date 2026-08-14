#ifndef SAVANNAH_PLATFORMS_OPENGL_SHADER
#define SAVANNAH_PLATFORMS_OPENGL_SHADER

#include <savannah/proto-core.h>
#include <external/glad/include/glad/glad.h>
#include <external/glfw/glfw3.h>

#include <external/glm/glm.hpp> // for vec and mat template types

#include <savannah/shader.h>

namespace Savannah
{
	class OpenGLShader : public Shader {
	public:
		OpenGLShader(const std::string& name, const std::string& vertexSrc, const std::string& fragmentSrc, const std::string& computeSrc);
		OpenGLShader(const std::string& name, const std::string& filepath);
		virtual ~OpenGLShader();

		virtual void Bind() const override;
		virtual void Unbind() const override;
		// virtual void UploadUniformBuffer(std::string name, T buffer) override;
		// virtual void UploadUniformBuffer() override

		void UploadUniformInt(std::string name, const int& value) const override;
		void UploadUniformIntArray(std::string name, int* values, int count) const override;
		void UploadUniformFloat(std::string name, const float& value) const override;
		void UploadUniformFloatArray(std::string name, float* values, int count) const override;

		void UploadUniformVec2(std::string name, const glm::vec2& vector) const override;
		void UploadUniformVec3(std::string name, const glm::vec3& vector) const override;
		void UploadUniformVec4(std::string name, const glm::vec4& vector) const override;

		void UploadUniformMat3(std::string name, const glm::mat3& matrix) const override;
		void UploadUniformMat4(std::string name, const glm::mat4& matrix) const override;


		//protected:
		//	uint32_t m_RendererId;
	private:
		std::unordered_map<GLenum, std::string> Preprocess(std::string& source);
		void Compile(const std::unordered_map<GLenum, std::string>& shaderSources);
	};
}

#endif
