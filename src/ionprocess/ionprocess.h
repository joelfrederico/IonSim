#ifndef __LIBIONPROCESS_H_INCLUDED__
#define __LIBIONPROCESS_H_INCLUDED__

#include <string>

std::string _getion(unsigned int step);
int index(int i, int j, int ncols);
int index(int i, int j, int k, int ncols, int nrows);
int makehist(std::string filename, int xbins, int step, unsigned long long *&hist, long &histsize, int &n_field_z);
int dealloc_hist(unsigned long long *&hist);
int get_n_field_z(std::string filename, int step);

#endif
