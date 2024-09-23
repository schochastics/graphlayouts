#define R_NO_REMAP
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>  // for NULL

#include <fstream>
#include <stdexcept>

#include "QuadCensus.h"

using namespace std;
using namespace oaqc;

static SEXP c_to_r(const unsigned int* const rMap, const unsigned long* orbit,
                   unsigned int size, unsigned long orbitCount) {
  const unsigned long len = size * orbitCount;
  SEXP rvec = PROTECT(Rf_allocVector(REALSXP, len));
  // set dimensions
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(dim)[0] = size;
  INTEGER(dim)[1] = orbitCount;
  Rf_setAttrib(rvec, R_DimSymbol, dim);
  // copy (cast and transpose) data
  double* data = REAL(rvec);
  unsigned int j = 0;

  for (unsigned int o = 0; o < orbitCount; ++o) {
    for (unsigned int i = 0; i < size; ++i) {
      unsigned int pos = i;
      if (rMap != NULL) {
        pos = rMap[i];
      }
      data[j++] = (double)orbit[pos * orbitCount + o];
    }
  }
  UNPROTECT(2);
  return rvec;
}

static void write_to_file(string fName, const unsigned int* const rMap,
                          const unsigned long* orbit, unsigned int size,
                          unsigned long orbitCount) {
  ofstream writer;
  // open the writer
  writer.open(fName.c_str());
  if (!writer.is_open()) {
    throw ios_base::failure("cannot open " + fName);
  }
  for (unsigned int o = 0; o < orbitCount - 1; ++o) {
    writer << "orbit_" << o << ";";
  }
  writer << "orbit_" << (orbitCount - 1) << std::endl;

  for (unsigned int i = 0; i < size; ++i) {
    unsigned int pos = i;
    if (rMap != NULL) {
      pos = rMap[i];
    }
    for (unsigned int o = 0; o < orbitCount - 1; ++o) {
      writer << orbit[pos * orbitCount + o] << ";";
    }
    writer << orbit[pos * orbitCount + orbitCount - 1] << std::endl;
  }
  // flush and close
  writer.flush();
  writer.close();
}

static void write_results(SEXP& a_value, SEXP& a_names, unsigned int& sIndex,
                          const unsigned int n, const unsigned int m,
                          QuadCensus& qc, string& filePrefix, string rType) {
  if (!filePrefix.empty()) {
    write_to_file(filePrefix + "_n_orbits_" + rType + ".csv", qc.getMapping(),
                  qc.nOrbits(), n, qc.getNOrbitCount());
    write_to_file(filePrefix + "_e_orbits_" + rType + ".csv", NULL,
                  qc.eOrbits(), m, qc.getEOrbitCount());
  }
  SET_STRING_ELT(a_names, sIndex, Rf_mkChar(("n_orbits_" + rType).c_str()));
  SET_VECTOR_ELT(a_value, sIndex,
                 c_to_r(qc.getMapping(), qc.nOrbits(), n, qc.getNOrbitCount()));
  ++sIndex;
  SET_STRING_ELT(a_names, sIndex, Rf_mkChar(("e_orbits_" + rType).c_str()));
  SET_VECTOR_ELT(a_value, sIndex,
                 c_to_r(NULL, qc.eOrbits(), m, qc.getEOrbitCount()));
  ++sIndex;
}

extern "C" SEXP entry(SEXP a_n, SEXP a_edges, SEXP a_freqFlag, SEXP a_file) {
  const unsigned int n = INTEGER(a_n)[0];
  const unsigned int m = Rf_length(a_edges) / 2;
  const int* edges = INTEGER(a_edges);
  string filePrefix(CHAR(STRING_ELT(a_file, 0)));
  const bool wNonIndFreq = LOGICAL(a_freqFlag)[0];

  unsigned int res_size;
  if (wNonIndFreq) {
    res_size = 4;
  } else {
    res_size = 2;
  }
  QuadCensus qc(n, m, edges);

  SEXP value = PROTECT(Rf_allocVector(VECSXP, res_size));  // VECSXP == list
  SEXP names = PROTECT(Rf_allocVector(STRSXP, res_size));

  unsigned int pos = 0;

  if (wNonIndFreq) {
    write_results(value, names, pos, n, m, qc, filePrefix, "non_ind");
  }

  qc.calcInducedFrequencies();

  write_results(value, names, pos, n, m, qc, filePrefix, "ind");

  Rf_setAttrib(value, R_NamesSymbol, names);

  UNPROTECT(2);
  return value;
}

#define CALLDEF(name, n) {#name, (DL_FUNC) & name, n}

static const R_CallMethodDef CallEntries[] = {CALLDEF(entry, 4),
                                              {NULL, NULL, 0}};

extern "C" void R_init_oaqc(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
