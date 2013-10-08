#include "typedefs.h"

#include <xmmintrin.h>
#include <emmintrin.h>

#include <assert.h>

void ConvertVideoFrame420ToRGB(
	const th_info *tinfo,
	const th_ycbcr_buffer ycbcr,
	unsigned char* pixels
)
{
	// some constant definitions 

	const float single_yoffset = 16.0f;
	const float single_yexcursion = 219.0f;
	const float single_cboffset = 128.0f;
	const float single_cbexcursion = 224.0f;
	const float single_croffset = 128.0f;
	const float single_crexcursion = 224.0f;

	const float kr = 0.299f;
	const float kb = 0.114f;

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

		const unsigned char* ydata = yplane.data;
		const unsigned char* cbdata = cbplane.data;
		const unsigned char* crdata = crplane.data;

		const int ystride = yplane.stride;
		const int cbstride = cbplane.stride;
		const int crstride = crplane.stride;

		const __m128 yoffset = _mm_set_ps1(-single_yoffset);
		const __m128 yexcursion = _mm_set_ps1(1.0f / single_yexcursion);
		const __m128 cboffset = _mm_set_ps1(-single_cboffset);
		const __m128 cbexcursion = _mm_set_ps1(1.0f / single_cbexcursion);
		const __m128 croffset = _mm_set_ps1(-single_croffset);
		const __m128 crexcursion = _mm_set_ps1(1.0f / single_crexcursion);

		const __m128 fr = _mm_set_ps1(255.0f * 2 * (1 - kr));
		const __m128 fb = _mm_set_ps1(255.0f * 2 * (1 - kb));

		const __m128 f1 = _mm_set_ps1(255.0f * (2 * (1 - kb) * kb / (1 - kb - kr)));
		const __m128 f2 = _mm_set_ps1(255.0f * (2 * (1 - kr) * kr / (1 - kb - kr)));

		const __m128 c255 = _mm_set_ps1(255.0f);

		for(int h = 0; h < height; ++h)
		{
			for(int w = 0; w < width; w += 16)
			{
				const __m128i yIn = _mm_loadu_si128((const __m128i*)(ydata + h*ystride + w));
				// assumption is that there is only one pixel in the cb/cr plane per 4 pixels (2x2) in the y plane
				const __m128i cbIn = _mm_loadu_si128((const __m128i*)(cbdata + h/2*cbstride + w/2));
				const __m128i crIn = _mm_loadu_si128((const __m128i*)(crdata + h/2*crstride + w/2));

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