#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_version.h>

void warn_handler (const char * reason, const char * file, int line, int gsl_errno) {
  warn("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
  return;
}

MODULE = Math::EasyGSL	PACKAGE = Math::EasyGSL

BOOT:
printf("Hello from the bootstrap!\n"); /* XXX-FIXME */
gsl_set_error_handler (&warn_handler);

void
GSL_VERSION()
	PPCODE:
		XPUSHs(sv_2mortal(newSVpv(GSL_VERSION,0)));

void
GSL_MAJOR_VERSION()
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(GSL_MAJOR_VERSION)));

void
GSL_MINOR_VERSION()
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(GSL_MINOR_VERSION)));
