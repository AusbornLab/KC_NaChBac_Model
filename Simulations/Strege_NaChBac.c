/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__na_bac_strege
#define _nrn_initial _nrn_initial__na_bac_strege
#define nrn_cur _nrn_cur__na_bac_strege
#define _nrn_current _nrn_current__na_bac_strege
#define nrn_jacob _nrn_jacob__na_bac_strege
#define nrn_state _nrn_state__na_bac_strege
#define _net_receive _net_receive__na_bac_strege 
#define states states__na_bac_strege 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define gmax_columnindex 0
#define ina _p[1]
#define ina_columnindex 1
#define ina2 _p[2]
#define ina2_columnindex 2
#define m _p[3]
#define m_columnindex 3
#define h _p[4]
#define h_columnindex 4
#define ena _p[5]
#define ena_columnindex 5
#define Dm _p[6]
#define Dm_columnindex 6
#define Dh _p[7]
#define Dh_columnindex 7
#define v _p[8]
#define v_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_tauh(void);
 static void _hoc_taum(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_na_bac_strege", _hoc_setdata,
 "tauh_na_bac_strege", _hoc_tauh,
 "taum_na_bac_strege", _hoc_taum,
 0, 0
};
#define tauh tauh_na_bac_strege
#define taum taum_na_bac_strege
 extern double tauh( _threadargsprotocomma_ double );
 extern double taum( _threadargsprotocomma_ double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_na_bac_strege", "S/cm2",
 "ina_na_bac_strege", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"na_bac_strege",
 "gmax_na_bac_strege",
 0,
 "ina_na_bac_strege",
 "ina2_na_bac_strege",
 0,
 "m_na_bac_strege",
 "h_na_bac_strege",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gmax = 0.12;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Strege_NaChBac_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 na_bac_strege Strege_NaChBac.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
    Dm = ( 1.0 / ( 1.0 + exp ( ( - 47.20712608 - v ) / 8.171491024 ) ) - m ) / taum ( _threadargscomma_ v ) ;
   Dh = ( 1.0 / ( 1.0 + exp ( ( - 56.68326979 - v ) / - 6.078198844 ) ) - h ) / tauh ( _threadargscomma_ v ) ;
    }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taum ( _threadargscomma_ v ) )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauh ( _threadargscomma_ v ) )) ;
   return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
     m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taum ( _threadargscomma_ v ))))*(- ( ( ( ( 1.0 ) / ( 1.0 + exp ( ( - 47.20712608 - v ) / 8.171491024 ) ) ) ) / taum ( _threadargscomma_ v ) ) / ( ( ( ( - 1.0 ) ) ) / taum ( _threadargscomma_ v ) ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauh ( _threadargscomma_ v ))))*(- ( ( ( ( 1.0 ) / ( 1.0 + exp ( ( - 56.68326979 - v ) / - 6.078198844 ) ) ) ) / tauh ( _threadargscomma_ v ) ) / ( ( ( ( - 1.0 ) ) ) / tauh ( _threadargscomma_ v ) ) - h) ;
    }
  return 0;
}
 
double taum ( _threadargsprotocomma_ double _lv ) {
   double _ltaum;
  _ltaum = - 125253.86 / ( 1.0 + exp ( - 0.06 * ( _lv - - 187.78 ) ) ) + 125258.35 ;
    
return _ltaum;
 }
 
static void _hoc_taum(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  taum ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double tauh ( _threadargsprotocomma_ double _lv ) {
   double _ltauh;
  _ltauh = 36.11944425 / ( 1.0 + exp ( ( - 24.58946004 - _lv ) / - 5.15842944 ) ) + 165.144429082 ;
    
return _ltauh;
 }
 
static void _hoc_tauh(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  tauh ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   m = 1.0 / ( 1.0 + exp ( ( - 47.20712608 - v ) / 8.171491024 ) ) ;
   h = 1.0 / ( 1.0 + exp ( ( - 56.68326979 - v ) / - 6.078198844 ) ) ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ina = gmax * pow( m , 3.0 ) * h * ( v - ena ) ;
   ina2 = pow( m , 3.0 ) * h ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "Strege_NaChBac.mod";
static const char* nmodl_file_text = 
  ": NaChBac channel model based on Strege et al., 2023.\n"
  ": Model developed by: Anthony Moreno-Sanchez \n"
  "\n"
  "\n"
  ": Set units\n"
  "UNITS {\n"
  "   (S)  = (siemens)\n"
  "   (mV) = (millivolt)\n"
  "   (mA) = (milliamp)\n"
  "}\n"
  "\n"
  ": Declare the channel\n"
  ": the suffix declaration in the below section refers to the name that this mod file is called by in python/hoc code (see line 20 of the python script)\n"
  "NEURON {\n"
  "   SUFFIX na_bac_strege\n"
  "   USEION na READ ena WRITE ina\n"
  "   RANGE gmax, ina, m , h, ina2\n"
  "}\n"
  "\n"
  ": Set default parameters\n"
  ": this is where you adjust gmax\n"
  "PARAMETER {\n"
  "   gmax = 0.12 (S/cm2) : the value to the left is the value in Neuron's original hh.mod file\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "   ena (mV)\n"
  "   ina (mA/cm2)\n"
  "   v  (mV)\n"
  "   ina2 ()\n"
  "}\n"
  "\n"
  "STATE { m h }\n"
  "\n"
  ": Calculate activation and current\n"
  ": this is where you adjust degrees for activation/inactivation\n"
  "BREAKPOINT {\n"
  "   SOLVE states METHOD cnexp\n"
  "   ina = gmax * m^3 * h * (v - ena)\n"
  "   ina2 =  m^3 * h \n"
  "}\n"
  "\n"
  ": Initialize activation to steady state\n"
  ": this is where steady state curves are defined\n"
  "INITIAL { \n"
  "\n"
  "   : values extracted from Strege et al. are tuned by cmaes\n"
  "   m = 1.0 / (1.0 + exp( (-47.20712608-v) / 8.171491024))\n"
  "   h = 1.0 / (1.0 + exp( (-56.68326979-v) / -6.078198844 ))\n"
  "}\n"
  "\n"
  ": Calculate change in activation \n"
  ": equations should have same values as ones above, make sure that the derivative states you are using matches the steady-state equations you are using \n"
  ": this means that the values in your m' should match the values in your m, and same with your h' and h\n"
  "DERIVATIVE states {\n"
  "   UNITSOFF\n"
  "\n"
  "   m' = (1.0 / (1.0 + exp( (-47.20712608-v) / 8.171491024 )) - m) / taum(v)\n"
  "   h' = (1.0 / (1.0 + exp( (-56.68326979-v) / -6.078198844 )) - h) / tauh(v)\n"
  "\n"
  "   UNITSON\n"
  "}\n"
  "\n"
  ": Function for activation time constant\n"
  ": taum is defined here\n"
  "FUNCTION taum(v (mV)) (/ms) {\n"
  "   UNITSOFF\n"
  "\n"
  "   taum = -125253.86 / (1 + exp(-0.06 * (v - -187.78))) + 125258.35\n"
  "\n"
  "\n"
  "   UNITSON\n"
  "}\n"
  "\n"
  ": Function for inactivation time constant\n"
  ": tauh is defined here\n"
  "FUNCTION tauh(v (mV)) (/ms) {\n"
  "   UNITSOFF\n"
  "   \n"
  "   tauh = 36.11944425 / ( 1 + exp((-24.58946004-v)/-5.15842944) ) + 165.144429082\n"
  "\n"
  "   UNITSON\n"
  "}\n"
  ;
#endif
