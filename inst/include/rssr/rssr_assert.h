#ifndef RSSRASSERT_H
#define RSSRASSERT_H
void pkassert (const char *msg, const char *file, int line);
#define rassert(EX)
//#define rassert(EX) (void)((EX) || (pkassert (#EX, __FILE__, __LINE__),0))

#endif
