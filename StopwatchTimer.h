#ifndef STOPWATCHTIME_H
#define STOPWATCHTIME_H
#include <windows.h>
#include <stdio.h>

#define  USE_MICROSECTIMER

#ifdef USE_MICROSECTIMER
struct MicroSecondTimer
{
	LARGE_INTEGER m_liPerfFreq;
	LARGE_INTEGER m_liPerfStart;
	LARGE_INTEGER m_liPerfEnd;
	SYSTEMTIME startT;
	SYSTEMTIME endT;
	SYSTEMTIME elastedT;
	std::vector<SYSTEMTIME> pauseT;
	std::vector<SYSTEMTIME> restartT;
	double m_liPerfElasted;
	int index;
	int fileTypeindex;
	//FILE * logFpsFile;
	std::vector<LARGE_INTEGER>  m_liPerfPause;
	std::vector<LARGE_INTEGER> m_liPerfRestart;
	void setid(int id)
	{
		index = id;
		fileTypeindex = 1;
		if (!QueryPerformanceFrequency(&m_liPerfFreq))
		{
			std::cout << "this PC doesn't support this timer" << std::endl;
			return;
		}
	}


	void startMyTimer()
	{
#ifdef USE_MILLITIMER
		GetLocalTime(&startT);
#endif


		QueryPerformanceCounter(&m_liPerfStart);
	}
	void endMyTimer()
	{
		QueryPerformanceCounter(&m_liPerfEnd);
		m_liPerfElasted = (((m_liPerfEnd.QuadPart -
			m_liPerfStart.QuadPart) * 1000000) / (double)m_liPerfFreq.QuadPart);

		for (unsigned int i = 0; i<m_liPerfRestart.size(); i++)
		{
			//elastedT.wSecond -=(restartT[i].wSecond - pauseT[i].wSecond);
			//int dist = restartT[i].wSecond - pauseT[i].wSecond;
			//if (dist>0)
			//	elastedT.wMilliseconds = elastedT.wMilliseconds - ((int)restartT[i].wMilliseconds + 1000 * dist - pauseT[i].wMilliseconds);
			//else 
			m_liPerfElasted = m_liPerfElasted - ((m_liPerfRestart[i].QuadPart - m_liPerfPause[i].QuadPart) * 1000000) / (double)m_liPerfFreq.QuadPart;

		}
		if (m_liPerfRestart.size()>0)
		{

			m_liPerfRestart.clear();
			m_liPerfPause.clear();
			std::vector<LARGE_INTEGER> swap1, swap2;
			swap1.swap(m_liPerfRestart);
			swap2.swap(m_liPerfPause);
		}

#ifdef USE_MILLITIMER
GetLocalTime(&endT);
//elastedT.wSecond = endT.wSecond - startT.wSecond;

elastedT.wMilliseconds = endT.wMilliseconds - startT.wMilliseconds;
if ((endT.wSecond - startT.wSecond)>0)
{
	int dist = endT.wSecond - startT.wSecond;
	elastedT.wMilliseconds = endT.wMilliseconds + 1000 * dist - startT.wMilliseconds;
}
if (restartT.size() > pauseT.size())
elastedT.wMilliseconds = pauseT[pauseT.size() - 1].wMilliseconds - startT.wMilliseconds;
for (int i = 0; i<restartT.size(); i++)
{
	//elastedT.wSecond -=(restartT[i].wSecond - pauseT[i].wSecond);
	int dist = restartT[i].wSecond - pauseT[i].wSecond;
	if (dist>0)
		elastedT.wMilliseconds = elastedT.wMilliseconds - ((int)restartT[i].wMilliseconds + 1000 * dist - pauseT[i].wMilliseconds);
	else elastedT.wMilliseconds = elastedT.wMilliseconds - (restartT[i].wMilliseconds - pauseT[i].wMilliseconds);

}
if (pauseT.size() > 0)
{
	restartT.clear();
	pauseT.clear();
	vector<SYSTEMTIME> swap1, swap2;
	swap1.swap(restartT);
	swap2.swap(pauseT);
}

#endif
//printTimer();
	}
	void pauseMyTimer()
	{
		LARGE_INTEGER tmp;
		QueryPerformanceCounter(&tmp);
		m_liPerfPause.push_back(tmp);
#ifdef USE_MILLITIMER
		SYSTEMTIME tmp1;
		GetLocalTime(&tmp1);
		pauseT.push_back(tmp1);
#endif
	}
	void restartMyTimer()
	{
		LARGE_INTEGER tmp;
		QueryPerformanceCounter(&tmp);
		m_liPerfRestart.push_back(tmp);
#ifdef USE_MILLITIMER
		SYSTEMTIME tmp1;
		GetLocalTime(&tmp1);
		restartT.push_back(tmp1);
#endif
	}
	void printTimer()
	{
		std::cout << "MicroTimer" << index << " went " << m_liPerfElasted << std::endl;
#ifdef USE_MILLITIMER
		printf("elasted:%02d,%03d\n", elastedT.wSecond, elastedT.wMilliseconds);
#endif
	}
	/*void printToFileElastedTime()
	{
		fprintf(logFpsFile, "%5.3f\n", m_liPerfElasted*0.001f);
	}*/
	float getDt()
	{
		return (float)m_liPerfElasted * 0.001f;
	}
	float getFPS()
	{
		if (m_liPerfElasted != 0)
		{
			float time = (1.0f / ((float)m_liPerfElasted * 0.001f * 0.001f));

			return   time;
		}
		else return 0.0f;
	}
};

static MicroSecondTimer mt, mtUpdate, mtDraw, mtTest;


#endif

#endif