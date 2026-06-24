#pragma once

#include <cstdio>
#include <filesystem>
#include <string>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#define NULL_DEVICE "NUL"
#define BENCH_DUP _dup
#define BENCH_DUP2 _dup2
#define BENCH_FILENO _fileno
#define BENCH_OPEN _open
#define BENCH_CLOSE _close
#define BENCH_WRONLY _O_WRONLY
#else
#include <fcntl.h>
#include <unistd.h>
#define NULL_DEVICE "/dev/null"
#define BENCH_DUP dup
#define BENCH_DUP2 dup2
#define BENCH_FILENO fileno
#define BENCH_OPEN open
#define BENCH_CLOSE close
#define BENCH_WRONLY O_WRONLY
#endif

namespace benchutil {

// Redirects stdout (and optionally stderr) to NULL_DEVICE for the lifetime of the
// object, restoring the original streams on destruction. Unlike the previous
// freopen("/dev/tty") approach, this saves/restores the real file descriptors with
// dup/dup2, so it works when there is no controlling terminal (CI, redirected runs).
class ScopedStdoutSuppressor {
public:
   explicit ScopedStdoutSuppressor(bool suppress_stderr = false)
   {
      std::fflush(stdout);
      m_savedStdout = BENCH_DUP(BENCH_FILENO(stdout));
      m_devnull = BENCH_OPEN(NULL_DEVICE, BENCH_WRONLY);
      if (m_devnull >= 0)
         BENCH_DUP2(m_devnull, BENCH_FILENO(stdout));
      if (suppress_stderr && m_devnull >= 0) {
         std::fflush(stderr);
         m_savedStderr = BENCH_DUP(BENCH_FILENO(stderr));
         BENCH_DUP2(m_devnull, BENCH_FILENO(stderr));
      }
   }

   ~ScopedStdoutSuppressor()
   {
      std::fflush(stdout);
      if (m_savedStdout >= 0) {
         BENCH_DUP2(m_savedStdout, BENCH_FILENO(stdout));
         BENCH_CLOSE(m_savedStdout);
      }
      if (m_savedStderr >= 0) {
         std::fflush(stderr);
         BENCH_DUP2(m_savedStderr, BENCH_FILENO(stderr));
         BENCH_CLOSE(m_savedStderr);
      }
      if (m_devnull >= 0)
         BENCH_CLOSE(m_devnull);
   }

   ScopedStdoutSuppressor(const ScopedStdoutSuppressor &) = delete;
   ScopedStdoutSuppressor &operator=(const ScopedStdoutSuppressor &) = delete;

private:
   int m_savedStdout = -1;
   int m_savedStderr = -1;
   int m_devnull = -1;
};

// Sum the sizes of all files in the current directory whose name contains `pattern`.
inline std::size_t GetTotalFileSize(const std::string &pattern)
{
   std::size_t total = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos)
         total += std::filesystem::file_size(entry.path());
   }
   return total;
}

// Remove all files in the current directory whose name contains `pattern`.
inline void CleanupFiles(const std::string &pattern)
{
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos)
         std::remove(entry.path().string().c_str());
   }
}

} // namespace benchutil
