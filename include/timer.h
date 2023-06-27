//
// Created by Zaiqiang Wu on 2023/06/18.
//

#ifndef A3_FINITE_ELEMENTS_3D_TIMER_H
#define A3_FINITE_ELEMENTS_3D_TIMER_H
#include<chrono>
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
class Timer
{
public:
    Timer()
    {
        elapsed_time=0.f;
    }
    void begin()
    {
        t_begin=high_resolution_clock::now();
    }
    float get_elapsed_time()
    {
        elapsed_time=(float)std::chrono::duration_cast<milliseconds>(high_resolution_clock::now()-t_begin).count();
        return elapsed_time;
    }

private:
    float elapsed_time;
    high_resolution_clock::time_point t_begin;
};

#endif //A3_FINITE_ELEMENTS_3D_TIMER_H
