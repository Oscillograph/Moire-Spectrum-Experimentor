#include <savannah/platforms/opengl/opengl_shader.h>
#include <savannah/utils/file-io.h> // FileIO::GetRawText()
#include <glm/gtc/type_ptr.hpp> // glm::value_ptr()
#include <savannah/utils/debugger.h> // SAVANNAH_CORE_ERROR

namespace Savannah
{
	static GLenum ShaderTypeFromString(const std::string& type) {
		if (type == "vertex"){
			return GL_VERTEX_SHADER;
		}
		if ((type == "fragment") || (type == "pixel")){
			return GL_FRAGMENT_SHADER;
		}
		if (type == "geometry"){
			return GL_GEOMETRY_SHADER;
		}
		if (type == "compute"){
			return GL_COMPUTE_SHADER;
		}
		SAVANNAH_CORE_ASSERT(false, "Unknown shader type: {0}", type);
		return GL_NONE;
	}

	OpenGLShader::OpenGLShader(const std::string& name, const std::string& filepath){
		std::string source = FileIO::GetRawText(filepath);

		if (source.size() > 1) {
			auto shaderSources = Preprocess(source);
			Compile(shaderSources);
			SetName(name);
		}

		/*
		 * std::string result;
		 * std::ifstream in(filepath, std::ios::in | std::ios::binary);
		 * if (in) { // read line by line
		 *	in.seekg(0, std::ios::end);
		 *	result.resize(in.tellg());
		 *	in.seekg(0, std::ios::beg);
		 *	in.read(&result[0], result.size());
		 *	in.close();
	} else {
		SAVANNAH_CORE_ERROR("Couldn't open file: {0}", filepath);
	}
	*/
	}

	OpenGLShader::OpenGLShader(const std::string& name, const std::string& vertexSrc, const std::string& fragmentSrc){
		std::unordered_map<GLenum, std::string> shaderSources;
		shaderSources[GL_VERTEX_SHADER] = vertexSrc;
		shaderSources[GL_FRAGMENT_SHADER] = fragmentSrc;
		Compile(shaderSources);
		SetName(name);
	}

	OpenGLShader::~OpenGLShader(){
		glDeleteProgram(GetRendererID());
	}

	void OpenGLShader::Bind() const {
		glUseProgram(GetRendererID());
	}

	void OpenGLShader::Unbind() const {
		glUseProgram(0);
	}

	void OpenGLShader::UploadUniformInt(std::string name, const int& value) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform1i(location, value);
	}

	void OpenGLShader::UploadUniformIntArray(std::string name, int* values, int count) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform1iv(location, count, values);
	}

	void OpenGLShader::UploadUniformFloat(std::string name, const float& value) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform1f (location, value);
	}

	void OpenGLShader::UploadUniformFloatArray(std::string name, float* values, int count) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform1fv(location, count, values);
	}

	void OpenGLShader::UploadUniformVec2(std::string name, const glm::vec2& vector) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform2fv (location, 2, glm::value_ptr(vector));
	}

	void OpenGLShader::UploadUniformVec3(std::string name, const glm::vec3& vector) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform3fv (location, 3, glm::value_ptr(vector));
	}

	void OpenGLShader::UploadUniformVec4(std::string name, const glm::vec4& vector) const {
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniform4fv (location, 4, glm::value_ptr(vector));
	}

	void OpenGLShader::UploadUniformMat3(std::string name, const glm::mat3& matrix) const {
		// glProgramUniformMatrix4fv
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniformMatrix3fv(location, 1, GL_FALSE, glm::value_ptr(matrix));
	}

	void OpenGLShader::UploadUniformMat4(std::string name, const glm::mat4& matrix) const {
		// glProgramUniformMatrix4fv
		GLint location = glGetUniformLocation(GetRendererID(), name.c_str());
		glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(matrix));
	}

	std::unordered_map<GLenum, std::string> OpenGLShader::Preprocess(std::string& source){
		std::unordered_map<GLenum, std::string> shaderSources;

		const char* typeToken = "#type";
		size_t typeTokenLength = strlen(typeToken);
		size_t pos = source.find(typeToken, 0);
		while (pos != std::string::npos) {
			size_t eol = source.find_first_of("\r\n", pos);
			SAVANNAH_CORE_ASSERT((eol != std::string::npos), "Syntax error");
			size_t begin = pos + typeTokenLength + 1; // only one space allowed
			std::string type = source.substr(begin, eol - begin);

			size_t nextLinePos = source.find_first_not_of("\r\n", eol);
			pos = source.find(typeToken, nextLinePos);
			shaderSources[ShaderTypeFromString(type)] = source.substr(
				nextLinePos,
				pos - (nextLinePos == std::string::npos ? (source.size() - 1) : nextLinePos)
			);
		}

		return shaderSources;
	}

	void OpenGLShader::Compile(const std::unordered_map<GLenum, std::string>& shaderSources){
		SetRendererID(glCreateProgram());
		GLuint program = GetRendererID();

		std::vector<GLuint> glShaderIDs;
		glShaderIDs.resize(shaderSources.size());
		//std::array<GLuint, 2> glShaderIDs;
		//int glShaderIDsIndex = 0;

		for (auto& key_value : shaderSources) {
			GLenum type = key_value.first;
			std::string GLSLsource = key_value.second;

			GLuint shader = glCreateShader(type);

			// Send the vertex shader source code to GL
			// Note that std::string's .c_str is NULL character terminated.
			const GLchar *source = (const GLchar *)GLSLsource.c_str();
			glShaderSource(shader, 1, &source, 0);

			// Compile the vertex shader
			glCompileShader(shader);

			GLint isCompiled = 0;
			glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
			if(isCompiled == GL_FALSE) {
				GLint maxLength = 0;
				glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

				// The maxLength includes the NULL character
				std::vector<GLchar> infoLog(maxLength);
				glGetShaderInfoLog(shader, maxLength, &maxLength, &infoLog[0]);

				// We don't need the shader anymore.
				glDeleteShader(shader);

				// Use the infoLog as you see fit.
				SAVANNAH_CORE_ERROR("{0}", infoLog.data());
				SAVANNAH_CORE_ASSERT(false, "Shader compilation failed.");
				break;
			}

			// Attach our shaders to our program
			glAttachShader(program, shader);
			glShaderIDs.push_back(shader);
			//glShaderIDs[glShaderIDsIndex] = shader;
			//glShaderIDsIndex++;
		}

		// Vertex and fragment shaders are successfully compiled.
		// Now time to link them together into a program.
		// Get a program object.
		m_RendererId = program;

		// Link our program
		glLinkProgram(program);

		// Note the different functions here: glGetProgram* instead of glGetShader*.
		GLint isLinked = 0;
		glGetProgramiv(program, GL_LINK_STATUS, (int *)&isLinked);
		if (isLinked == GL_FALSE) {
			GLint maxLength = 0;
			glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

			// The maxLength includes the NULL character
			std::vector<GLchar> infoLog(maxLength);
			glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);

			// We don't need the program anymore.
			glDeleteProgram(program);
			// Don't leak shaders either.
			for (auto shaderID : glShaderIDs){
				glDeleteShader(shaderID);
			}

			// Use the infoLog as you see fit.
			SAVANNAH_CORE_ERROR("{0}", infoLog.data());
			SAVANNAH_CORE_ASSERT(false, "OpenGL program compilation failed due to linking error.");

			// In this simple program, we'll just leave
			return;
		}

		// Always detach shaders after a successful link.
		for (auto shaderID : glShaderIDs){
			glDetachShader(program, shaderID);
		}

		//if (shaderSources.size() <= 2){
		//} else {
		//	SAVANNAH_CORE_ERROR("BSE supports only 2 shaders for now");
		//}
		//std::vector<GLuint> glShaderIDs;
		//glShaderIDs.reserve(shaderSources.size());
	}
}
