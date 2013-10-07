#include <fstream>
#include <iostream>
#include <iomanip> 
#include <Windows.h>

using namespace std;

typedef bool bool_t;
typedef unsigned char ubyte_t;
typedef double double_t;

typedef struct{
  /**The width of this plane.*/
  int            width;
  /**The height of this plane.*/
  int            height;
  /**The offset in bytes between successive rows.*/
  int            stride;
  /**A pointer to the beginning of the first row.*/
  unsigned char *data;
}th_img_plane;

typedef unsigned __int32 ogg_uint32_t;

typedef enum{
  TH_CS_UNSPECIFIED,
  TH_CS_ITU_REC_470M,
  TH_CS_ITU_REC_470BG,
  TH_CS_NSPACES
}th_colorspace;

typedef enum{
  TH_PF_420,
  TH_PF_RSVD,
  TH_PF_422,
  TH_PF_444,
  TH_PF_NFORMATS
}th_pixel_fmt;

typedef struct{
  unsigned char version_major;
  unsigned char version_minor;
  unsigned char version_subminor;
  ogg_uint32_t  frame_width;
  ogg_uint32_t  frame_height;
  ogg_uint32_t  pic_width;
  ogg_uint32_t  pic_height;
  ogg_uint32_t  pic_x;
  ogg_uint32_t  pic_y;
  ogg_uint32_t  fps_numerator;
  ogg_uint32_t  fps_denominator;
  ogg_uint32_t  aspect_numerator;
  ogg_uint32_t  aspect_denominator;
  th_colorspace colorspace;
  th_pixel_fmt  pixel_fmt;
  int           target_bitrate;
  int           quality;
  int           keyframe_granule_shift;
}th_info;

typedef th_img_plane th_ycbcr_buffer[3];

void read_function(th_info* ti, th_ycbcr_buffer yuv)
{
	ifstream file;
	file.open("test.dmp", ios_base::binary);

	file.read((char*)ti, sizeof(th_info));

	for (int i = 0; i < 3; ++i)
	{
		file.read((char*)&yuv[i], sizeof(th_img_plane));
		size_t data_length = yuv[i].height * yuv[i].stride;
		yuv[i].data = new unsigned char[data_length];
		file.read((char*)yuv[i].data, data_length);

		if (!file)
		{
			cout << "error:" << endl;
		}
	}

	file.close();
}

bool WriteBmp(
	const char* path,
	const int width,
	const int height,
	const unsigned char* image
)
{
	ofstream file;

	file.open(path, ios_base::binary);
	if (!file)
	{
		return false;
	}

	BITMAPFILEHEADER bmfh;
	BITMAPINFOHEADER info;
	memset ( &bmfh, 0, sizeof (BITMAPFILEHEADER ) );
	memset ( &info, 0, sizeof (BITMAPINFOHEADER ) );

	bmfh.bfType = 0x4d42; // 0x4d42 = 'BM'
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	bmfh.bfOffBits = 0x36;

	info.biSize = sizeof(BITMAPINFOHEADER);
	info.biWidth = width;
	info.biHeight = height;
	info.biPlanes = 1;
	info.biBitCount = 24;
	info.biCompression = BI_RGB;
	info.biSizeImage = 0;
	info.biXPelsPerMeter = 0x0ec4;
	info.biYPelsPerMeter = 0x0ec4;
	info.biClrUsed = 0;	
	info.biClrImportant = 0; 

	if (!file.write(reinterpret_cast<const char*>(&bmfh), sizeof(BITMAPFILEHEADER )))
	{	
		return false;
	}

	if (!file.write(reinterpret_cast<const char*>(&info), sizeof(BITMAPINFOHEADER )))
	{	
		return false;
	}

	if (!file.write(reinterpret_cast<const char*>(image), width*height*3))
	{	
		return false;
	}

	return true;
}

class Profiler
{
public:
	Profiler();
	inline ~Profiler();

	inline void Start();
	inline void Stop();
	inline void Reset();
	inline const double_t GetTime() const;

protected:

	bool_t is_running_;
	double_t time_;

	LONGLONG start_;
	LONGLONG running_time_;
};

inline Profiler::Profiler()
:	is_running_(false),
	time_(0.0),
	start_(0),
	running_time_(0)
{
}

inline Profiler::~Profiler()
{
}

inline void Profiler::Start()
{
	QueryPerformanceCounter((LARGE_INTEGER *)&start_);

	is_running_ = true;
}

inline void Profiler::Stop()
{
	LONGLONG Stop;
	QueryPerformanceCounter((LARGE_INTEGER *)&Stop);
	running_time_ += (Stop - start_);

	is_running_ = false;
}

inline void Profiler::Reset()
{
	LONGLONG freq;
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	double_t double_running_time = (double_t)running_time_;
	double_t double_freq = (double_t)freq;
	time_ = double_running_time / double_freq;
	running_time_ = 0;
}

inline const double_t Profiler::GetTime() const
{
	return time_;
}

#include <xmmintrin.h>
#include <emmintrin.h>

#include <tmmintrin.h>

#include <assert.h>

void ConvertVideoFrame420ToRGB_(
	const th_info *tinfo,
	const th_ycbcr_buffer ycbcr,
	unsigned char* pixels
)
{
	const float single_yoffset = 16.0f;
	const float single_yexcursion = 219.0f;
	const float single_cboffset = 128.0f;
	const float single_cbexcursion = 224.0f;
	const float single_croffset = 128.0f;
	const float single_crexcursion = 224.0f;
	const float kr = 0.299f;
	const float kb = 0.114f;

	__m128 yoffset = _mm_set_ps1(single_yoffset);
	__m128 yexcursion = _mm_set_ps1(1.0f / single_yexcursion);
	__m128 c255 = _mm_set_ps1(255.0f);
	__m128 c0 = _mm_set_ps1(0.0f);

	const int w = tinfo->pic_width;
	const int half_w = w/2;
	const int h = tinfo->pic_height;
	if (pixels)
	{
		unsigned char *dst = pixels;
		const int ystride = ycbcr[0].stride;
		const int cbstride = ycbcr[1].stride;
		const int crstride = ycbcr[2].stride;
		const int dststride = 3*w;
		const int yoff = (tinfo->pic_x & ~1) + ystride * (tinfo->pic_y & ~1);
		const int cboff = (tinfo->pic_x / 2) + (cbstride) * (tinfo->pic_y / 2);
		const unsigned char *py = ycbcr[0].data + yoff;
		const unsigned char *pcb = ycbcr[1].data + cboff;
		const unsigned char *pcr = ycbcr[2].data + cboff;
		int posx, posy;
		for (posy = 0; posy < h; posy += 2)
		{
			const unsigned char* py_1 = py;
			const unsigned char* py_2 = py + ystride;

			unsigned char* dst_1 = dst + (h - posy - 1) * dststride;
			unsigned char* dst_2 = dst + (h - posy - 2) * dststride;

			for (posx = 0; posx < half_w; ++posx)
			{
				// http://www.theora.org/doc/Theora.pdf, 1.1 spec,
				// chapter 4.2 (Y'CbCr -> Y'PbPr -> R'G'B')
				// These constants apparently work for NTSC _and_ PAL/SECAM.

				float yf[4] = {
					static_cast<float>(py_1[2*posx+0]),
					static_cast<float>(py_1[2*posx+1]),
					static_cast<float>(py_2[2*posx+0]),
					static_cast<float>(py_2[2*posx+1]),
				};

				__m128 y = _mm_load_ps(yf);

				y = _mm_sub_ps(y, yoffset);
				y = _mm_mul_ps(y, yexcursion);

				const float pb = (((float) pcb[posx]) - single_cboffset) / single_cbexcursion;
				const float pr = (((float) pcr[posx]) - single_croffset) / single_crexcursion;

				// Convert YPbPr to RGB
				__m128 tmp_r = _mm_set_ps1(2.0f * (1.0f - kr) * pr);
				__m128 tmp_g = _mm_set_ps1((2.0f * (((1.0f - kb) * kb) / (1.0f - kb - kr))) * pb + (2.0f * (((1.0f - kr) * kr) / (1.0f - kb - kr))) * pr);
				__m128 tmp_b = _mm_set_ps1(2.0f * (1.0f - kb) * pb);

				__m128 r = _mm_mul_ps(_mm_add_ps(y, tmp_r), c255);
				__m128 g = _mm_mul_ps(_mm_sub_ps(y, tmp_g), c255);
				__m128 b = _mm_mul_ps(_mm_add_ps(y, tmp_b), c255);

				// clamp to range [0, 255]
				r = _mm_max_ps(c0, _mm_min_ps(c255, r));
				g = _mm_max_ps(c0, _mm_min_ps(c255, g));
				b = _mm_max_ps(c0, _mm_min_ps(c255, b));

				float rf[4];
				float gf[4];
				float bf[4];

				_mm_store_ps(rf, r);
				_mm_store_ps(gf, g);
				_mm_store_ps(bf, b);

				*(dst_1++) = static_cast<ubyte_t>(bf[0]);
				*(dst_1++) = static_cast<ubyte_t>(gf[0]);
				*(dst_1++) = static_cast<ubyte_t>(rf[0]);

				*(dst_1++) = static_cast<ubyte_t>(bf[1]);
				*(dst_1++) = static_cast<ubyte_t>(gf[1]);
				*(dst_1++) = static_cast<ubyte_t>(rf[1]);

				*(dst_2++) = static_cast<ubyte_t>(bf[2]);
				*(dst_2++) = static_cast<ubyte_t>(gf[2]);
				*(dst_2++) = static_cast<ubyte_t>(rf[2]);

				*(dst_2++) = static_cast<ubyte_t>(bf[3]);
				*(dst_2++) = static_cast<ubyte_t>(gf[3]);
				*(dst_2++) = static_cast<ubyte_t>(rf[3]);
			} // for

			// adjust to the start of the next line.
			py += 2*ystride;
			pcb += cbstride;
			pcr += crstride;
			//dst += 2*dststride;
		} // for
	} // if
}

void ConvertVideoFrame420ToRGB(
	const th_info *tinfo,
	const th_ycbcr_buffer ycbcr,
	unsigned char* pixels
)
{
	if (pixels)
	{
		const th_img_plane yplane = ycbcr[0];
		const th_img_plane cbplane = ycbcr[1];
		const th_img_plane crplane = ycbcr[2];

		const int width = tinfo->pic_width;
		const int height = tinfo->pic_height;
		const int wh = width*height;
		
		assert(wh == yplane.width*yplane.height);
		assert(width % 16 == 0);
		assert(cbplane.width * 2 == yplane.width);
		assert(crplane.width * 2 == yplane.width);

		const float single_yoffset = 16.0f;
		const float single_yexcursion = 219.0f;
		const float single_cboffset = 128.0f;
		const float single_cbexcursion = 224.0f;
		const float single_croffset = 128.0f;
		const float single_crexcursion = 224.0f;
		const float single_kr = 0.299f;
		const float single_kb = 0.114f;

		const __m128 yoffset = _mm_set_ps1(-single_yoffset);
		const __m128 yexcursion = _mm_set_ps1(1.0f / single_yexcursion);
		const __m128 cboffset = _mm_set_ps1(-single_cboffset);
		const __m128 cbexcursion = _mm_set_ps1(1.0f / single_cbexcursion);
		const __m128 croffset = _mm_set_ps1(-single_croffset);
		const __m128 crexcursion = _mm_set_ps1(1.0f / single_crexcursion);

		const __m128 kr = _mm_set_ps1(single_kr);
		const __m128 kb = _mm_set_ps1(single_kb);
		const __m128 c255 = _mm_set_ps1(255.0f);
		const __m128 c1 = _mm_set_ps1(1.0f);
		const __m128 c2 = _mm_set_ps1(2.0f);

		const __m128 fr = _mm_mul_ps(c255, _mm_mul_ps(c2, _mm_sub_ps(c1, kr)));
		const __m128 fb = _mm_mul_ps(c255, _mm_mul_ps(c2, _mm_sub_ps(c1, kb)));

		const __m128 f1 = _mm_mul_ps(c255, _mm_mul_ps(c2, _mm_div_ps(
				_mm_mul_ps(_mm_sub_ps(c1, kb), kb),
				_mm_sub_ps(_mm_sub_ps(c1, kb), kr))));
		const __m128 f2 = _mm_mul_ps(c255, _mm_mul_ps(c2, _mm_div_ps(
				_mm_mul_ps(_mm_sub_ps(c1, kr), kr),
				_mm_sub_ps(_mm_sub_ps(c1, kb), kr))));

		for(int h = 0; h < height; ++h)
		{
			for(int w = 0; w < width; w += 16)
			{
				const __m128i yIn = _mm_loadu_si128((const __m128i*)(yplane.data + h*width + w));
				// assumption is that there is only one pixel in the cb/cr plane per 4 pixels (2x2) in the y plane
				const __m128i cbIn = _mm_loadu_si128((const __m128i*)(cbplane.data + h/2*cbplane.stride + w/2));
				const __m128i crIn = _mm_loadu_si128((const __m128i*)(crplane.data + h/2*crplane.stride + w/2));

				// yIn ep8 -> ps

				const __m128i yInlo = _mm_unpacklo_epi8((yIn), _mm_setzero_si128());
				const __m128i yInHi = _mm_unpackhi_epi8((yIn), _mm_setzero_si128());
				const __m128i yIn1 = _mm_unpacklo_epi16(yInlo, _mm_setzero_si128());
				const __m128i yIn4 = _mm_unpackhi_epi16(yInlo, _mm_setzero_si128());
				const __m128i yIn8 = _mm_unpacklo_epi16(yInHi, _mm_setzero_si128());
				const __m128i yIn12 = _mm_unpackhi_epi16(yInHi, _mm_setzero_si128());

				const __m128 yIn1ps = _mm_cvtepi32_ps(yIn1);
				const __m128 yIn2ps = _mm_cvtepi32_ps(yIn4);
				const __m128 yIn3ps = _mm_cvtepi32_ps(yIn8);
				const __m128 yIn4ps = _mm_cvtepi32_ps(yIn12);

				// cbIn ep8 -> ps

				const __m128i cbInExp = _mm_unpacklo_epi8(cbIn, cbIn);
				const __m128i cbInlo = _mm_unpacklo_epi8(cbInExp, _mm_setzero_si128());
				const __m128i cbInHi = _mm_unpackhi_epi8(cbInExp, _mm_setzero_si128());
				const __m128i cbIn1 = _mm_unpacklo_epi16(cbInlo, _mm_setzero_si128());
				const __m128i cbIn4 = _mm_unpackhi_epi16(cbInlo, _mm_setzero_si128());
				const __m128i cbIn8 = _mm_unpacklo_epi16(cbInHi, _mm_setzero_si128());
				const __m128i cbIn12 = _mm_unpackhi_epi16(cbInHi, _mm_setzero_si128());

				const __m128 cbIn1ps = _mm_cvtepi32_ps(cbIn1);
				const __m128 cbIn2ps = _mm_cvtepi32_ps(cbIn4);
				const __m128 cbIn3ps = _mm_cvtepi32_ps(cbIn8);
				const __m128 cbIn4ps = _mm_cvtepi32_ps(cbIn12);

				// crIn ep8 -> ps

				const __m128i crInExp = _mm_unpacklo_epi8(crIn, crIn);
				const __m128i crInlo = _mm_unpacklo_epi8(crInExp, _mm_setzero_si128());
				const __m128i crInHi = _mm_unpackhi_epi8(crInExp, _mm_setzero_si128());
				const __m128i crIn1 = _mm_unpacklo_epi16(crInlo, _mm_setzero_si128());
				const __m128i crIn4 = _mm_unpackhi_epi16(crInlo, _mm_setzero_si128());
				const __m128i crIn8 = _mm_unpacklo_epi16(crInHi, _mm_setzero_si128());
				const __m128i crIn12 = _mm_unpackhi_epi16(crInHi, _mm_setzero_si128());

				const __m128 crIn1ps = _mm_cvtepi32_ps(crIn1);
				const __m128 crIn2ps = _mm_cvtepi32_ps(crIn4);
				const __m128 crIn3ps = _mm_cvtepi32_ps(crIn8);
				const __m128 crIn4ps = _mm_cvtepi32_ps(crIn12);

				// map [0..255] to [-1/2..+1/2] resp. [0..1]
			
				const __m128 yOut1ps = _mm_mul_ps(_mm_add_ps(yIn1ps, yoffset), yexcursion);
				const __m128 yOut2ps = _mm_mul_ps(_mm_add_ps(yIn2ps, yoffset), yexcursion);
				const __m128 yOut3ps = _mm_mul_ps(_mm_add_ps(yIn3ps, yoffset), yexcursion);
				const __m128 yOut4ps = _mm_mul_ps(_mm_add_ps(yIn4ps, yoffset), yexcursion);

				const __m128 cbOut1ps = _mm_mul_ps(_mm_add_ps(cbIn1ps, cboffset), cbexcursion);
				const __m128 cbOut2ps = _mm_mul_ps(_mm_add_ps(cbIn2ps, cboffset), cbexcursion);
				const __m128 cbOut3ps = _mm_mul_ps(_mm_add_ps(cbIn3ps, cboffset), cbexcursion);
				const __m128 cbOut4ps = _mm_mul_ps(_mm_add_ps(cbIn4ps, cboffset), cbexcursion);

				const __m128 crOut1ps = _mm_mul_ps(_mm_add_ps(crIn1ps, croffset), crexcursion);
				const __m128 crOut2ps = _mm_mul_ps(_mm_add_ps(crIn2ps, croffset), crexcursion);
				const __m128 crOut3ps = _mm_mul_ps(_mm_add_ps(crIn3ps, croffset), crexcursion);
				const __m128 crOut4ps = _mm_mul_ps(_mm_add_ps(crIn4ps, croffset), crexcursion);

				// do the actual conversion math (on range 0..255/-127..127 instead or 0..1/-1/2..+1/2

				const __m128 y1_255 = _mm_mul_ps(c255, yOut1ps);
				const __m128 y2_255 = _mm_mul_ps(c255, yOut2ps);
				const __m128 y3_255 = _mm_mul_ps(c255, yOut3ps);
				const __m128 y4_255 = _mm_mul_ps(c255, yOut4ps);

				const __m128 r1_1 = _mm_add_ps(y1_255, _mm_mul_ps(fr, crOut1ps));
				const __m128 r2_1 = _mm_add_ps(y2_255, _mm_mul_ps(fr, crOut2ps));
				const __m128 r3_1 = _mm_add_ps(y3_255, _mm_mul_ps(fr, crOut3ps));
				const __m128 r4_1 = _mm_add_ps(y4_255, _mm_mul_ps(fr, crOut4ps));

				const __m128 g1_1 = _mm_sub_ps(_mm_sub_ps(y1_255, _mm_mul_ps(f1, cbOut1ps)), _mm_mul_ps(f2, crOut1ps));
				const __m128 g2_1 = _mm_sub_ps(_mm_sub_ps(y2_255, _mm_mul_ps(f1, cbOut2ps)), _mm_mul_ps(f2, crOut2ps));
				const __m128 g3_1 = _mm_sub_ps(_mm_sub_ps(y3_255, _mm_mul_ps(f1, cbOut3ps)), _mm_mul_ps(f2, crOut3ps));
				const __m128 g4_1 = _mm_sub_ps(_mm_sub_ps(y4_255, _mm_mul_ps(f1, cbOut4ps)), _mm_mul_ps(f2, crOut4ps));

				const __m128 b1_1 = _mm_add_ps(y1_255, _mm_mul_ps(fb, cbOut1ps));
				const __m128 b2_1 = _mm_add_ps(y2_255, _mm_mul_ps(fb, cbOut2ps));
				const __m128 b3_1 = _mm_add_ps(y3_255, _mm_mul_ps(fb, cbOut3ps));
				const __m128 b4_1 = _mm_add_ps(y4_255, _mm_mul_ps(fb, cbOut4ps));

				// clip to 255
				const __m128 r1 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, r1_1));
				const __m128 r2 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, r2_1));
				const __m128 r3 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, r3_1));
				const __m128 r4 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, r4_1));

				const __m128 g1 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, g1_1));
				const __m128 g2 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, g2_1));
				const __m128 g3 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, g3_1));
				const __m128 g4 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, g4_1));

				const __m128 b1 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, b1_1));
				const __m128 b2 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, b2_1));
				const __m128 b3 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, b3_1));
				const __m128 b4 = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(c255, b4_1));

				// multiplex rgb channels

#define rgb_multiplex(no) \
				const __m128 rgb##no##_1 = _mm_shuffle_ps( \
					_mm_shuffle_ps(b##no##, g##no##, _MM_SHUFFLE(0, 0, 0, 0)),  \
					_mm_shuffle_ps(r##no##, b##no##, _MM_SHUFFLE(1, 1, 0, 0)),  \
					_MM_SHUFFLE(2, 0, 2, 0)); \
				const __m128 rgb##no##_2 = _mm_shuffle_ps( \
					_mm_shuffle_ps(g##no##, r##no##, _MM_SHUFFLE(1, 1, 1, 1)),  \
					_mm_shuffle_ps(b##no##, g##no##, _MM_SHUFFLE(2, 2, 2, 2)),  \
					_MM_SHUFFLE(2, 0, 2, 0)); \
				const __m128 rgb##no##_3 = _mm_shuffle_ps( \
					_mm_shuffle_ps(r##no##, b##no##, _MM_SHUFFLE(3, 3, 2, 2)),  \
					_mm_shuffle_ps(g##no##, r##no##, _MM_SHUFFLE(3, 3, 3, 3)),  \
					_MM_SHUFFLE(2, 0, 2, 0)); 

				rgb_multiplex(1);
				rgb_multiplex(2);
				rgb_multiplex(3);
				rgb_multiplex(4);

#undef rgb_multiplex

				// pack 32bit -> 8bit

				const __m128i pack1l = _mm_packs_epi32(_mm_cvtps_epi32(rgb1_1), _mm_cvtps_epi32(rgb1_2));
				const __m128i pack1h = _mm_packs_epi32(_mm_cvtps_epi32(rgb1_3), _mm_cvtps_epi32(rgb2_1));
				const __m128i pack1 = _mm_packus_epi16(pack1l, pack1h);

				const __m128i pack2l = _mm_packs_epi32(_mm_cvtps_epi32(rgb2_2), _mm_cvtps_epi32(rgb2_3));
				const __m128i pack2h = _mm_packs_epi32(_mm_cvtps_epi32(rgb3_1), _mm_cvtps_epi32(rgb3_2));
				const __m128i pack2 = _mm_packus_epi16(pack2l, pack2h);

				const __m128i pack3l = _mm_packs_epi32(_mm_cvtps_epi32(rgb3_3), _mm_cvtps_epi32(rgb4_1));
				const __m128i pack3h = _mm_packs_epi32(_mm_cvtps_epi32(rgb4_2), _mm_cvtps_epi32(rgb4_3));
				const __m128i pack3 = _mm_packus_epi16(pack3l, pack3h);

				// and finally store in output

				_mm_storeu_si128((__m128i*)(pixels + ((wh-width)*3) - h*width*3 + w*3 + 0*16), pack1);
				_mm_storeu_si128((__m128i*)(pixels + ((wh-width)*3) - h*width*3 + w*3 + 1*16), pack2);
				_mm_storeu_si128((__m128i*)(pixels + ((wh-width)*3) - h*width*3 + w*3 + 2*16), pack3);
			}
		}

	} // if
}

int main(const int argc, const char* argv[])
{
	th_info ti;
	th_ycbcr_buffer yuv;
	read_function(&ti, yuv);

	unsigned char* buffer = new unsigned char[3*ti.pic_width*ti.pic_height];

	for(int j = 0; j < 5; ++j)
	{
		Profiler p;
		p.Start();
		for(int i = 0; i < 2000; ++i)
		{
			ConvertVideoFrame420ToRGB(&ti, yuv, buffer);
		}
		p.Stop();
		p.Reset();

		cout << "Time taken: " << p.GetTime() << endl;
	}

	WriteBmp("out.bmp", ti.pic_width, ti.pic_height, buffer);

	return 0;
}
