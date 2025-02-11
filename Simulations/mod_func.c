#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Strege_NaChBac_reg();
extern void _ks_gunay_reg();
extern void _nap_reg();
extern void _nat_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Strege_NaChBac.mod");
fprintf(stderr," ks_gunay.mod");
fprintf(stderr," nap.mod");
fprintf(stderr," nat.mod");
fprintf(stderr, "\n");
    }
_Strege_NaChBac_reg();
_ks_gunay_reg();
_nap_reg();
_nat_reg();
}
