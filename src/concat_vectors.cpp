#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List concat_vectors(List list1, List list2) {
  int n = list1.size();
  List result(n);
  
  for(int i = 0; i < n; i++) {
    CharacterVector v1 = list1[i];
    CharacterVector v2 = list2[i];
    
    // Concatenate v2 to v1 and store in result
    CharacterVector combined(v1.size() + v2.size());
    std::copy(v1.begin(), v1.end(), combined.begin());
    std::copy(v2.begin(), v2.end(), combined.begin() + v1.size());

    result[i] = combined;
  }
  
  return result;
}