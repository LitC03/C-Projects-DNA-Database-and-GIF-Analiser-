// tictoc.h
#ifndef TICTOC_H
#define TICTOC_H

#include <iostream>
#include <chrono>

using namespace std;

class TicToc
{
private:
  chrono::steady_clock::time_point time_start;
  chrono::steady_clock::time_point time_end;
  chrono::duration<double> time_elapsed_seconds;

public:
  TicToc();

  void tic(); // start timer
  void toc(); // end timer

  friend ostream& operator<<(ostream& os, const TicToc& time_in);
};

#endif // TICTOC_H