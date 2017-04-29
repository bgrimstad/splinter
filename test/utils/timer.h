/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TIMER_H
#define SPLINTER_TIMER_H

#include "chrono"

namespace SPLINTER {

class Timer
{
public:
    Timer()
            : t_start(std::chrono::steady_clock::now()),
              t_duration(std::chrono::duration<double>(0)),
              running(false)
    {
    }

    void start()
    {
        if (!running)
        {
            t_start = std::chrono::steady_clock::now();
            running = true;
        }
    }

    void stop()
    {
        if (running)
        {
            auto tNow = std::chrono::steady_clock::now();
            auto dur = (tNow-t_start);
            t_duration += dur;
            running = false;
        }
    }

    void reset()
    {
        t_duration = std::chrono::duration<double>(0);
        running = false;
    }

    long get_seconds()
    {
        return std::chrono::duration_cast<std::chrono::seconds> (getDuration()).count();
    }

    long get_milli_seconds()
    {
        return std::chrono::duration_cast<std::chrono::milliseconds> (getDuration()).count();
    }

    long get_micro_seconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds> (getDuration()).count();
    }

private:
    std::chrono::steady_clock::time_point t_start;
    std::chrono::duration<double> t_duration;

    bool running;

    std::chrono::duration<double> getDuration()
    {
        if (running)
        {
            auto tNow = std::chrono::steady_clock::now();
            auto dur = (tNow-t_start);
            return (t_duration + dur);
        }
        else
        {
            return t_duration;
        }
    }

};

} // namespace SPLINTER

#endif //SPLINTER_TIMER_H
