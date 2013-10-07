#include <fstream>
#include <iostream>
#include <iomanip> 
#include <Windows.h>

#include "typedefs.h"

using namespace std;

typedef bool bool_t;
typedef unsigned char ubyte_t;
typedef double double_t;

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

void ConvertVideoFrame420ToRGB(
	const th_info *tinfo,
	const th_ycbcr_buffer ycbcr,
	unsigned char* pixels
);

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
