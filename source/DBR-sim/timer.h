// Copied from https://gist.github.com/mcleary/b0bf4fa88830ff7c882d
// Created by mcleary
#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <cmath>

/**
 * @class Timer
 * @brief High-resolution timer for measuring elapsed time.
 * 
 * Provides methods to start, stop, and query elapsed time in milliseconds or seconds.
 * Uses std::chrono for precise timing measurements.
 */
class Timer
{
public:
    /**
     * @brief Start the timer.
     * 
     * Records the current time as the start time and sets the running state to true.
     */
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }

    /**
     * @brief Stop the timer.
     * 
     * Records the current time as the end time and sets the running state to false.
     */
    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }

    /**
     * @brief Get elapsed time in milliseconds.
     * 
     * @return double Elapsed time in milliseconds since start() was called.
     *         If timer is running, returns time elapsed so far.
     *         If timer is stopped, returns time between start() and stop().
     */
    double elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;

        if (m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }

        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }

    /**
     * @brief Get elapsed time in seconds.
     * 
     * @return double Elapsed time in seconds since start() was called.
     */
    double elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};