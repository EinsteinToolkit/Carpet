#ifndef PTHREAD_WRAPPER_HPP
#define PTHREAD_WRAPPER_HPP

#include <cctk.h>

#include <functional>
#include <string>
#include <utility>

#if HAVE_CAPABILITY_PTHREADS
#include <pthread.h>
#endif

namespace CarpetSimulationIO {

#if HAVE_CAPABILITY_PTHREADS

class pthread_wrapper_t {
  typedef std::function<void()> fun_t;

  std::string name;
  fun_t fun;
  pthread_t thread;

  static void *run(void *arg_);
  void *run1();

public:
  pthread_wrapper_t() = delete;
  pthread_wrapper_t(const pthread_wrapper_t &) = delete;
  pthread_wrapper_t(pthread_wrapper_t &&) = delete;
  pthread_wrapper_t &operator=(const pthread_wrapper_t &) = delete;
  pthread_wrapper_t &operator=(pthread_wrapper_t &&) = delete;

  pthread_wrapper_t(const std::string &name, const fun_t &fun);
  ~pthread_wrapper_t();
};

#else

struct pthread_wrapper_t {
  pthread_wrapper_t(const std::string &name, const fun_t &fun) { fun(); }
};

#endif
}

#endif // #ifndef PTHREAD_WRAPPER_HPP
