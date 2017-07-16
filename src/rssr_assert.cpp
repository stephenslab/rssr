#include "rssr.h"

void pkassert (const char *msg, const char *file, int line) {
  char buffer [200];
  snprintf( buffer, 200, "Assert: %s at %s line #%d\n", msg, file,
	    line );
  ::Rf_error( buffer );
}
