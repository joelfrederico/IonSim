#ifndef __LIBIONPROCESS_H_INCLUDED__
#define __LIBIONPROCESS_H_INCLUDED__

#include <string>
#include <vector>

std::string _getion(unsigned int step);
int index(int i, int j, int ncols);
int index(int i, int j, int k, int ncols, int nrows);
int makehist(std::string filename, int xbins, int step, std::vector<unsigned long long> &hist, long &histsize, int &n_field_z);

#endif
