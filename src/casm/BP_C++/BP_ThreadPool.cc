#ifndef BP_ThreadPool_CC
#define BP_ThreadPool_CC

#include <stdio.h>
#include <cstdlib>
#include "casm/BP_C++/BP_ThreadPool.hh"

namespace BP {
  //public:

  // Construct a ThreadPool with 'Nthreads'
  //
  ThreadPool::ThreadPool(int Nthreads) {
    // create thread array
    N = Nthreads;
    thread = new pthread_t[N];

    // init the mutex and cond
    pthread_mutex_init(&queue_mutex, NULL);
    pthread_mutex_init(&work_count_mutex, NULL);
    pthread_cond_init(&work_avail_cond, NULL);
    pthread_cond_init(&queue_empty_cond, NULL);
    pthread_cond_init(&work_finished_cond, NULL);

    // init the work_count
    pthread_mutex_lock(&work_count_mutex);
    work_count = 0;
    pthread_mutex_unlock(&work_count_mutex);

    // create the worker threads
    for(int i = 0; i < N; i++) {
      pthread_create(&thread[i], NULL, worker_thread, (void *) this);
    }

    // the queue begins unlocked, but empty
  }

  // Destructor, waits until all work is finished before completing
  //
  ThreadPool::~ThreadPool() {
    int rc;
    void *status;

    // send the KILL signal to all workers
    //   (they will continue working until all tasks are complete)
    //
    for(int i = 0; i < N; i++) {
      pthread_mutex_lock(&queue_mutex);
      queue.push_back(Work(NULL, NULL, KILL));
      pthread_cond_signal(&work_avail_cond);
      pthread_mutex_unlock(&queue_mutex);

    }


    // wait for all workers to finish
    for(int i = 0; i < N; i++) {
      rc = pthread_join(thread[i], &status);
      if(rc) {
        printf("ERROR; return code from thread %d is %d\n", i, rc);
        std::exit(-1);
      }
    }

    // delete thread
    delete [] thread;

    // delete the mutex and cond
    pthread_mutex_destroy(&queue_mutex);
    pthread_mutex_destroy(&work_count_mutex);
    pthread_cond_destroy(&work_avail_cond);
    pthread_cond_destroy(&queue_empty_cond);
    pthread_cond_destroy(&work_finished_cond);
  };

  // Add work to do. Work is performed in a first-in first-out manner
  //
  //   Arguments:
  //     void* (*fp)(void *): function pointer to function that will execute
  //     void* farg: pointer to argument given to function pointer
  //
  //   ** Caution: arguments are passed as pointers, so you must make  **
  //   **          certain they remain valid until all work is done!   **
  //
  void ThreadPool::add_work(void *(*fp)(void *), void *farg) {
    // lock the queue
    pthread_mutex_lock(&queue_mutex);

    // add a task
    queue.push_back(Work(fp, farg, WORK));

    // signal that work is available
    pthread_cond_signal(&work_avail_cond);

    // unlock the queue
    pthread_mutex_unlock(&queue_mutex);

  };

  // Block until all work is finished
  //   (and do not allow adding more work until current work is finished)
  //
  void ThreadPool::finish() {
    // lock the queue
    pthread_mutex_lock(&queue_mutex);

    while(!queue.empty()) {
      pthread_cond_wait(&queue_empty_cond, &queue_mutex);
    }

    // lock the work_count
    pthread_mutex_lock(&work_count_mutex);

    while(work_count > 0) {
      pthread_cond_wait(&work_finished_cond, &work_count_mutex);
    }

    // unlock the work_count
    pthread_mutex_unlock(&work_count_mutex);

    // unlock the queue
    pthread_mutex_unlock(&queue_mutex);
  }

  // Check if finished
  bool ThreadPool::is_finished() {
    bool result = false;

    // lock the queue
    pthread_mutex_lock(&queue_mutex);

    // lock the work_count
    pthread_mutex_lock(&work_count_mutex);

    if(work_count == 0 && queue.empty())
      result = true;

    // unlock the work_count
    pthread_mutex_unlock(&work_count_mutex);

    // unlock the queue
    pthread_mutex_unlock(&queue_mutex);

    return result;
  };

  //private:

  void ThreadPool::worker() {
    do {
      // collect a task from the threadpool queue
      Work task = request_work();

      // if given work, execute the task
      if(task.is_work()) {
        task.execute();
        finish_work();
      }

      // if given KILL signal, exit do-while loop
      if(task.is_kill())
        break;
    }
    while(true);
  }

  // static Function run in individual threads
  void *ThreadPool::worker_thread(void *arg) {
    ((ThreadPool *) arg)->worker();
    return NULL;
  }

  // Called by 'worker' to request work
  Work ThreadPool::request_work() {
    // lock the queue
    pthread_mutex_lock(&queue_mutex);

    // if no work, wait for it to come along
    while(queue.empty()) {
      pthread_cond_wait(&work_avail_cond, &queue_mutex);
    }

    Work task = queue.front();
    queue.pop_front();

    if(queue.empty())
      pthread_cond_signal(&queue_empty_cond);

    // lock the work_count
    pthread_mutex_lock(&work_count_mutex);
    work_count++;

    // unlock the work_count
    pthread_mutex_unlock(&work_count_mutex);


    // unlock the queue
    pthread_mutex_unlock(&queue_mutex);

    return task;
  }

  // Called by 'worker' when task is complete
  void ThreadPool::finish_work() {
    // lock the work_count
    pthread_mutex_lock(&work_count_mutex);

    work_count--;

    if(work_count == 0)
      pthread_cond_signal(&work_finished_cond);

    // unlock the work_count
    pthread_mutex_unlock(&work_count_mutex);

  }

}

#endif // BP_ThreadPool_CC

