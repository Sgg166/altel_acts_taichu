#pragma once

#include <cstdio>
#include <csignal>

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>
#include <regex>
#include <iostream>
#include <future>

#include <unistd.h>

#include "myrapidjson.h"

namespace altel{
  class TelViewer{
  public:
    ~TelViewer();

    void pushObject(std::shared_ptr<JsonValue> dp);
    void popFrontObject();
    void clearObjects();
    std::shared_ptr<JsonValue>& frontObject();

    void stopAsyncLoop();
    void startAsyncLoop();
    uint64_t asyncLoop();

  private:
    std::vector<std::shared_ptr<JsonValue>> m_vec_ring_ev;
    std::shared_ptr<JsonValue> m_ring_end; // ring end is nullptr,therefore real data as nullptr should not go into the ring.

    uint64_t m_size_ring{200000};
    std::atomic_uint64_t m_count_ring_write;
    std::atomic_uint64_t m_count_ring_read;
    std::atomic_uint64_t m_hot_p_read;
    bool m_is_async_thread_running{false};
    std::future<uint64_t> m_fut_async_thread;
  };
}
