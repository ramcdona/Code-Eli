/*********************************************************************************
* Copyright (c) 2026 OpenVSP
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
********************************************************************************/

#ifndef eli_util_scratch_hpp
#define eli_util_scratch_hpp

#include <deque>
#include <cstddef>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace util
  {
    // Reusable per-thread scratch storage.
    //
    // Many hot routines need a temporary working buffer (typically a small dynamically-sized Eigen
    // matrix) that they fill, use, and discard on every call.  Allocating that buffer per call is
    // pure churn.  This facility hands out a reusable buffer instead:
    //
    //     eli::util::scratch< Eigen::Matrix<T, R, C> > buf;   // borrow
    //     *buf = something;                                    // resize is a no-op once warm
    //     ... use *buf / buf-> ...                             // buffer is returned when buf dies
    //
    // A single pool serves every routine that asks for a given buffer type, so a buffer grown by
    // one routine is reused by another -- the same scratch reused across different applications.
    // Each borrow returns a *distinct* buffer, so nested or recursive borrows are safe; buffers are
    // reused across borrows that do not overlap in time (the pool is a simple LIFO high-water
    // arena).  Because the borrowed matrix keeps its capacity between borrows, an Eigen resize back
    // to the same shape does not reallocate, so warm hot paths allocate nothing.
    //
    // The pool is thread_local, so concurrent evaluation from multiple threads is race-free without
    // any locking.

    template <typename Mat_>
    class scratch_pool
    {
      public:
        // The pool for the calling thread.
        static scratch_pool & instance()
        {
          static thread_local scratch_pool pool;
          return pool;
        }

        // Borrow the next free buffer, growing the arena (with stable addresses) if needed.
        Mat_ * acquire()
        {
          while ( m_buffers.size() <= m_top )
          {
            m_buffers.emplace_back();
          }
          return &m_buffers[ m_top++ ];
        }

        // Return the most recently borrowed buffer (strictly LIFO).
        void release()
        {
          --m_top;
        }

      private:
        std::deque< Mat_ > m_buffers;   // deque: growth never invalidates element addresses
        std::size_t m_top = 0;          // high-water cursor: index of the next free buffer
    };

    // RAII handle to a borrowed scratch buffer.  Behaves like a pointer to Mat_.
    template <typename Mat_>
    class scratch
    {
      public:
        scratch() : m_mat( scratch_pool<Mat_>::instance().acquire() ) {}
        ~scratch() { scratch_pool<Mat_>::instance().release(); }

        Mat_ & operator*()  const { return *m_mat; }
        Mat_ * operator->() const { return  m_mat; }
        Mat_ * get()        const { return  m_mat; }

        scratch( const scratch & ) = delete;
        scratch & operator=( const scratch & ) = delete;

      private:
        Mat_ * m_mat;
    };
  }
}

#endif
