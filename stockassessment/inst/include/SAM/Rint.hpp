
// Handle NA_integer on arm64
// Modified from: https://github.com/kaskr/adcomp/tree/naint

HEADER(
struct Rint {
  int i;
  inline Rint () { }
  inline Rint (double x) : i(IntegerFromReal(x)) {}
  inline operator int() const { return i; }
  // From src/main/coerce.c
  inline int IntegerFromReal(double x) {
    using std::isnan;
    if (ISNAN(x))
      return NA_INTEGER;
    else if (x >= INT_MAX+1. || x <= INT_MIN ) {
      return NA_INTEGER;
    }
    return (int) x;
  }
};
       )


// HEADER(
//        template<>
//        inline vector<int> asVector(SEXP x);
//        )

HEADER(
template<>
inline vector<int> asVector(SEXP x) {
  return asVector<Rint>(x).cast<int>();
}
       )

// HEADER(
//       template<>
//       inline matrix<int> asMatrix(SEXP x)
//        )

HEADER(
template<>
inline matrix<int> asMatrix(SEXP x) {
   return asMatrix<Rint>(x).cast<int>();
}
       )

// HEADER(
//        namespace tmbutils {
// 	 template<>
// 	 inline array<int> asArray(SEXP x);
//        }
//        )

HEADER(
namespace tmbutils {
template<>
inline array<int> asArray(SEXP x) {
  array<Rint> tmp = asArray<Rint>(x);
  return array<int>(tmp.vectorcopy.cast<int>(), tmp.dim);
}
}
       )
