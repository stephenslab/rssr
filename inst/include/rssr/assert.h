#ifndef DEBUG
#define ASSERT(x)
#else
#define ASSERT(x)                                 
if (! (x))                                        
{                                                 
  cout << "ERROR!! Assert " << #x << " failed\n"; 
  cout << " on line " << __LINE__  << "\n";       
  cout << " in file " << __FILE__ << "\n";        
}
#endif