#include <savannah/platforms/opengl/opengl_texture.h>

//#ifdef SAVANNAH_PLATFORM_WINDOWS
//	// for stb_image.h : Unicode filenames on Windows
//	#define STBI_WINDOWS_UTF8
//#endif
//#define STB_IMAGE_IMPLEMENTATION
#include <external/stb_image.h>

namespace Savannah
{
	OpenGLTexture2D::OpenGLTexture2D(const std::string& filename, bool savePixels)
		: m_Path(filename)
	{
		if (filename.size() > 0)
		{
			LoadTextureFromFile(m_Path, savePixels);
		}
	}
	
	OpenGLTexture2D::OpenGLTexture2D(unsigned char* pixels, int width, int height, int channels)
	{
		LoadTextureFromPixels(pixels, width, height, channels);
	}
	
	OpenGLTexture2D::~OpenGLTexture2D()
	{
		if (m_TextureID)
			glDeleteTextures(1, &m_TextureID);
		if (m_pixels != nullptr)
		{
			stbi_image_free(m_pixels);
		}
	}
	
	bool OpenGLTexture2D::LoadTextureFromFile(const std::string& filename, bool savePixels)
	{
		// Load from file
		bool result = false;
		int imageWidth, imageHeight, channels;
//		stbi_set_flip_vertically_on_load(1);
		stbi_uc* data = stbi_load(filename.c_str(), &imageWidth, &imageHeight, &channels, 0);
		
		if (data)
		{
			uint32_t internalFormat = GL_RGB8;
			uint32_t dataFormat = GL_RGB8;
			
			if (channels == 4){
				internalFormat = GL_RGBA8;
				dataFormat = GL_RGBA;
			} else if (channels == 3){
				internalFormat = GL_RGB8;
				dataFormat = GL_RGB;
			}
			
			m_InternalFormat = internalFormat;
			m_DataFormat = dataFormat;
			
			if ((internalFormat == 0) && ( dataFormat == 0))
			{
				printf("Image format is not supported\n");
			} else {
				m_TextureID = 0;
				m_ImageSize = {imageWidth, imageHeight};
				
				// Create a OpenGL texture identifier
				glGenTextures(1, &m_TextureID);
				glBindTexture(GL_TEXTURE_2D, m_TextureID);
				
				// What to do if the texture doesn't fit the place
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
				
				if (channels == 3)
				{
					glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
				}
				
				glTexImage2D(GL_TEXTURE_2D, 0, m_InternalFormat, imageWidth, imageHeight, 0, m_DataFormat, GL_UNSIGNED_BYTE, data);
				if (!savePixels)
				{
					stbi_image_free(data);
				} else {
					m_pixels = data;
				}
				result = true;
			}
		} else {
			printf("Couldn't load a texture from \"%s\"\n", filename.c_str());
		}
		
		return result;
	}
	
	bool OpenGLTexture2D::LoadTextureFromPixels(unsigned char* pixels, int width, int height, int channels)
	{
		bool result = false;
		
		if (pixels == nullptr)
		{
			printf("OpenGLTexture::LoadTextureFromPixels: no pixels received.\n");
			return result;
		}
		
		if ((width == 0) || (height == 0) || (channels == 0))
		{
			printf("OpenGLTexture::LoadTextureFromPixels: width, height and channels must have non-zero values.\n");
			return result;
		}
		
		if ((width > GL_MAX_TEXTURE_SIZE) || (height > GL_MAX_TEXTURE_SIZE))
		{
			printf("OpenGLTexture::LoadTextureFromPixels: width or height must not be greater than %d\n", GL_MAX_TEXTURE_SIZE);
			return result;
		}
		
		if (pixels)
		{
//            printf("Pixels address: ");
//            printf("%X\n", pixels);
			uint32_t internalFormat = GL_RGB8;
			uint32_t dataFormat = GL_RGB8;
			
			if (channels == 4){
				internalFormat = GL_RGBA8;
				dataFormat = GL_RGBA;
			} else if (channels == 3){
				internalFormat = GL_RGB8;
				dataFormat = GL_RGB;
			} else {
				printf("OpenGLTexture::LoadTextureFromPixels: only 4 or 3 channels supported.\n");
				return result;
			}
			
			m_InternalFormat = internalFormat;
			m_DataFormat = dataFormat;
			
			if ((internalFormat == 0) && ( dataFormat == 0))
			{
				printf("Image format is not supported\n");
			} else {
				m_TextureID = 0;
				m_ImageSize = {width, height};
				
				// Create a OpenGL texture identifier
				glGenTextures(1, &m_TextureID);
				glBindTexture(GL_TEXTURE_2D, m_TextureID);
				
				// What to do if the texture doesn't fit the place
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

				if (channels == 3)
				{
					glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
				}
			
//                printf("OpenGLTexture::LoadTextureFromPixels: Given pixel data in first 5x5 square:\n");
//                for (int k = 0; k < 5; ++k)
//                {
//                    for (int h = 0; h < 5; ++h)
//                    {
//                        printf("%X ", ((uint32_t*)pixels)[k*width + h]);
//                    }
//                    printf("\n");
//                }
                
				glTexImage2D(GL_TEXTURE_2D, 0, m_InternalFormat, width, height, 0, m_DataFormat, GL_UNSIGNED_BYTE, pixels);

				result = true;
//                printf("...done\n");
			}
		}
		
		return result;
	}
	
	GLuint OpenGLTexture2D::GetID()
	{
		return m_TextureID;
	}
	
	iVec2 OpenGLTexture2D::GetSize()
	{
		return m_ImageSize;
	}
	
	const iVec2 OpenGLTexture2D::GetSize() const
	{
		return m_ImageSize;
	}
	
	void OpenGLTexture2D::SetData(void* data, uint32_t size) {
		glTexImage2D(GL_TEXTURE_2D, 0, m_InternalFormat, m_ImageSize.x, m_ImageSize.y, 0, m_DataFormat, GL_UNSIGNED_BYTE, data);
		//BSE_CORE_INFO("glTexImage2D successful");
	}
	
	void OpenGLTexture2D::Bind(uint32_t slot) const {
		// bind to slot 0
		// glBindTextureUnit(slot, m_RendererID);
		// glBindTexture(slot, m_RendererID);
		glActiveTexture(GL_TEXTURE0 + slot); // for batch texture renderer
		// glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, m_TextureID);
		// BSE_INFO("Texture is bound to: {0}", GL_TEXTURE0 + slot);
	}
	
	void OpenGLTexture2D::Unbind(uint32_t slot) const {
		glActiveTexture(GL_TEXTURE0 + slot); // for batch texture renderer
		// glDisable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
		// BSE_INFO("Texture is unbound from: {0}", GL_TEXTURE0 + slot);
	}
}
