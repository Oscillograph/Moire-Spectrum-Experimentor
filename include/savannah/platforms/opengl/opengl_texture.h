#include <savannah/proto-core.h>
#include <external/glad/include/glad/glad.h>
#include <external/glfw/glfw3.h>

namespace Savannah
{
	class OpenGLTexture2D
	{
	public:
		OpenGLTexture2D(const std::string& filename = "", bool savePixels = false);
		OpenGLTexture2D(unsigned char* pixels = nullptr, int width = 0, int height = 0, int channels = 0);
		~OpenGLTexture2D();
		bool LoadTextureFromFile(const std::string& filename, bool savePixels = false);
		bool LoadTextureFromPixels(unsigned char* pixels = nullptr, int width = 0, int height = 0, int channels = 0);
		GLuint GetID();
		iVec2 GetSize();
		const iVec2 GetSize() const;
		void SetData(void* data, uint32_t size);
		void Bind(uint32_t slot) const;
		void Unbind(uint32_t slot) const;
		
		inline bool operator==(const OpenGLTexture2D* other) const
		{
			return (m_TextureID == ((OpenGLTexture2D*)other)->GetID());
		}
	private:
		std::string m_Path = "";
		uint32_t m_TextureID = 0;
		iVec2 m_ImageSize = {0, 0};
		unsigned char* m_pixels = nullptr;
		uint32_t m_InternalFormat = GL_RGB8;
		uint32_t m_DataFormat = GL_RGB8;
	};
}
