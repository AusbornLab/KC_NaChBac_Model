#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ks_gunay_reg(void);
extern void _nap_reg(void);
extern void _nat_reg(void);
extern void _Strege_NaChBac_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"ks_gunay.mod\"");
    fprintf(stderr, " \"nap.mod\"");
    fprintf(stderr, " \"nat.mod\"");
    fprintf(stderr, " \"Strege_NaChBac.mod\"");
    fprintf(stderr, "\n");
  }
  _ks_gunay_reg();
  _nap_reg();
  _nat_reg();
  _Strege_NaChBac_reg();
}

#if defined(__cplusplus)
}
#endif
