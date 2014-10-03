#ifndef OPERATORS_HH
#define OPERATORS_HH

// Transport (i.e., prolongation and restriction) operator types

enum operator_type
{
  op_error,                     // illegal operator type
  op_none,                      // do not transport
  op_sync,                      // transport only on the same level
                                // (error if called between levels)
  op_restrict,                  // restrict only, do not prolongate
  op_copy,                      // use simple copying for prolongation
                                // (needs only one time level)
  op_Lagrange,                  // Lagrange interpolation (standard)
  op_ENO,                       // use ENO stencils (for hydro)
  op_ENOVOL,                    // use ENO stencils, conserving mass (for hydro)
  op_WENO,                      // use WENO stencils (for hydro)
  op_TVD,                       // use TVD stencils (for hydro)
  op_Lagrange_monotone          // monotone Lagrange interpolation (for hydro)
};

#endif // OPERATORS_HH
