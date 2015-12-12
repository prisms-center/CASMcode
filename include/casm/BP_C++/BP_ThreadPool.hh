#ifndef BP_ThreadPool_HH
#define BP_ThreadPool_HH

#include <list>
#include <pthread.h>


// Simple multithreading
//   To link, include: -lpthread

namespace BP {
  enum WORKTAG {WORK, KILL};

  class Work {

    // function pointer that contains the 'work'
    void *(*fp)(void *);

    // argument to the function pointer
    void *farg;

    // tag saying if this is 'WORK' or 'KILL' signal
    WORKTAG tag;

  public:

    Work() {}

    Work(void *(*_fp)(void *), void *_farg, WORKTAG _tag): fp(_fp), farg(_farg), tag(_tag) {}

    // for now, no returning results
    void execute() {
      (*fp)(farg);
    }

    bool is_kill() {
      return (tag == KILL);
    }

    bool is_work() {
      return (tag == WORK);
    }
  };


  class ThreadPool {
    int N;
    pthread_t *thread;

    pthread_mutex_t queue_mutex;
    std::list< Work > queue;

    pthread_mutex_t work_count_mutex;
    int work_count;

    pthread_cond_t work_avail_cond;
    pthread_cond_t queue_empty_cond;
    pthread_cond_t work_finished_cond;


  public:

    // Construct a ThreadPool with 'Nthreads'
    //
    ThreadPool(int Nthreads);

    // Destructor, waits until all work is finished before completing
    //
    ~ThreadPool();

    // Add work to do. Work is performed in a first-in first-out manner
    //
    //   Arguments:
    //     void* (*fp)(void *): function pointer to function that will execute
    //     void* farg: pointer to argument given to function pointer
    //
    //   ** Caution: arguments are passed as pointers, so you must make  **
    //   **          certain they remain valid until all work is done!   **
    //
    void add_work(void *(*fp)(void *), void *farg);

    // Block until all work is finished
    //   (and do not allow adding more work until current work is finished)
    //
    void finish();

    // Check if finished
    bool is_finished();

  private:

    // Function run in individual threads
    void worker();

    // static Function that calls 'worker' for each thread
    static void *worker_thread(void *arg);

    // Called by 'worker' to request work
    Work request_work();

    // Called by 'worker' when task is complete
    void finish_work();

  };
}

#endif // BP_ThreadPool_HH

