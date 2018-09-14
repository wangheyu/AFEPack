/**
 * @file   Thread.cpp
 * @author Ruo Li
 * @date   Mon Feb 14 13:30:24 2005
 * 
 * @brief  for convenient using of POSIX thread in C++
 * 
 * 
 */

#ifdef MULTITHREAD

int n_thread = 1;

void setThread(unsigned int n)
{
  if (n < 1)
    n = 1;
  n_thread = n; 
};

unsigned int getThread()
{
  return n_thread;
};

#endif

/**
 * end of file
 * 
 */
