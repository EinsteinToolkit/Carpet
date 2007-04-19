#ifndef OPERATORS_HH
#define OPERATORS_HH

// Transport (i.e., prolongation and restriction) operator types

enum operator_type
{
  op_error,                     // illegal operator type
  op_none,                      // do not transport
  op_sync,                      // transport only on the same level
                                // (error if called between levels)
  op_copy,                      // use simple copying for prolongation
                                // (needs only one time level)
  op_Lagrange,                  // Lagrange interpolation (standard)
  op_ENO,                       // use ENO stencils (for hydro)
  op_WENO                       // use WENO stencils (for hydro)
};

#endif // OPERATORS_HH
