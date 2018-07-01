#pragma once

#ifdef __APPLE__ 
typedef __SIZE_TYPE__ size_t;
#endif

#include <string>
#include <stdio.h>  /* defines FILENAME_MAX */

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <direct.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

#if defined(_WIN32)
	class Lock {
		HANDLE hMutex;
	public:
		Lock()  { hMutex = CreateMutex (NULL, FALSE, NULL); }
		~Lock() { CloseHandle(hMutex); }
		void acquire() { WaitForSingleObject( hMutex, INFINITE ); }
		void release() { ReleaseMutex( hMutex ); }
	};
	static void yield(bool andSleep) { if (andSleep) Sleep(1); else SwitchToThread(); }
#else
	class Lock {
		pthread_mutex_t _mutex;
	public:
		Lock()  { _mutex = PTHREAD_MUTEX_INITIALIZER; }
		~Lock() {  }
		void acquire() { pthread_mutex_lock(&_mutex); }
		void release() { pthread_mutex_unlock(&_mutex); }
	};
	static void yield(bool andSleep) { if (andSleep) sleep(1); else sched_yield(); }
#endif



static std::string currentDirectory() {
	char buffer[FILENAME_MAX];
	if (!
#if defined(_WIN32)
		_getcwd
#else
		getcwd
#endif
			(buffer, FILENAME_MAX)) {
		return "<currentDirectory() failed>";
	}
	return buffer;
}

class ScopedAcquire{
public:
	Lock & lock;
	ScopedAcquire(Lock & lock) : lock(lock) { lock.acquire(); }
	~ScopedAcquire() { lock.release(); }
};
