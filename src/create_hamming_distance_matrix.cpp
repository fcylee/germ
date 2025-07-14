#include <Rcpp.h>
using namespace Rcpp;

// Helper function to calculate Hamming distance between two strings
int hamming_distance(const std::string& s1, const std::string& s2) {
  int dist = 0;
  for(size_t i = 0; i < s1.length(); ++i) {
    if(s1[i] != s2[i]) dist++;
  }
  return dist;
}

// Helper function to generate all k-mers
std::vector<std::string> generate_kmers(int k_len) {
  std::vector<std::string> kmers;
  std::vector<char> bases = {'A', 'C', 'G', 'T'};
  
  int total_kmers = std::pow(4, k_len);
  kmers.reserve(total_kmers + 1);
  
  // Mimic R's expand.grid behavior
  for(int i = 0; i < total_kmers; ++i) {
    std::string kmer(k_len, 'A');
    int temp = i;
    
    // Generate in expand.grid order (leftmost position cycles fastest)
    for(int pos = 0; pos < k_len; ++pos) {
      kmer[pos] = bases[temp % 4];
      temp /= 4;
    }
    kmers.push_back(kmer);
  }
  
  // Add N k-mer
  std::string n_kmer(k_len, 'N');
  kmers.push_back(n_kmer);
  
  return kmers;
}

//' Create scoring matrix based on hamming distance for kmers
//'
//' @param k_len an integer specifying the length of the k-mer
//' @param lambda an integer specifying the lambda for the exponential decay function. Set to NULL to use custom scale function
//' @param unweighted logical, if TRUE returns matrix of all 1s (except N vs N = 0)
//' @param scale_fun R function for scaling hamming distances (default: function(x) 1/(1+x^3))
//' @return NumericMatrix of k-mer similarity scores where element (i,j) represents similarity between k-mer i and k-mer j
//' @export
// [[Rcpp::export]]
NumericMatrix create_hamming_distance_matrix_cpp(int k_len, 
                                                Nullable<double> lambda = 1.0, 
                                                bool unweighted = false, 
                                                Nullable<Function> scale_fun = R_NilValue) {
  
  // Generate all k-mers
  std::vector<std::string> kmers = generate_kmers(k_len);
  int n_kmers = kmers.size();
  
  // Add row/column names first
  CharacterVector kmer_names(n_kmers);
  for(int i = 0; i < n_kmers; ++i) {
    kmer_names[i] = kmers[i];
  }
  
  // Handle unweighted case
  if(unweighted) {
    NumericMatrix result(n_kmers, n_kmers);
    
    // Fill with 1s
    for(int i = 0; i < n_kmers; ++i) {
      for(int j = 0; j < n_kmers; ++j) {
        result(i, j) = 1.0;
      }
    }
    
    // Set N k-mer vs N k-mer to 0
    result(n_kmers - 1, n_kmers - 1) = 0.0;
    
    rownames(result) = kmer_names;
    colnames(result) = kmer_names;
    return result;
  }
  
  // Calculate all pairwise Hamming distances
  NumericMatrix hamming_matrix(n_kmers, n_kmers);
  for(int i = 0; i < n_kmers; ++i) {
    for(int j = 0; j < n_kmers; ++j) {
      hamming_matrix(i, j) = hamming_distance(kmers[i], kmers[j]);
    }
  }
  
  NumericMatrix scaled_hamming_matrix(n_kmers, n_kmers);
  
  // Apply scaling based on priority: lambda > scale_fun > raw distances
  if(!Rf_isNull(lambda)) {
    // Use exponential decay function: exp(-lambda * x)
    double lambda_val = as<double>(lambda);
    for(int i = 0; i < n_kmers; ++i) {
      for(int j = 0; j < n_kmers; ++j) {
        double hamming_dist = hamming_matrix(i, j);
        scaled_hamming_matrix(i, j) = std::exp(-lambda_val * hamming_dist);
      }
    }
  } else if(!Rf_isNull(scale_fun)) {
    // Use custom scaling function
    Function scale_function = as<Function>(scale_fun);
    for(int i = 0; i < n_kmers; ++i) {
      for(int j = 0; j < n_kmers; ++j) {
        double hamming_dist = hamming_matrix(i, j);
        scaled_hamming_matrix(i, j) = as<double>(scale_function(hamming_dist));
      }
    }
  } else {
    // Return raw distances (no scaling)
    scaled_hamming_matrix = clone(hamming_matrix);
  }
  
  // Normalize to 0-1 range only if min != max
  double min_val = scaled_hamming_matrix(0, 0);
  double max_val = scaled_hamming_matrix(0, 0);
  
  for(int i = 0; i < n_kmers; ++i) {
    for(int j = 0; j < n_kmers; ++j) {
      if(scaled_hamming_matrix(i, j) < min_val) min_val = scaled_hamming_matrix(i, j);
      if(scaled_hamming_matrix(i, j) > max_val) max_val = scaled_hamming_matrix(i, j);
    }
  }
  
  if(min_val == max_val) {
    // All values are the same, set to 1
    for(int i = 0; i < n_kmers; ++i) {
      for(int j = 0; j < n_kmers; ++j) {
        scaled_hamming_matrix(i, j) = 1.0;
      }
    }
  } else {
    // Normalize: subtract min, then divide by max
    for(int i = 0; i < n_kmers; ++i) {
      for(int j = 0; j < n_kmers; ++j) {
        scaled_hamming_matrix(i, j) = (scaled_hamming_matrix(i, j) - min_val) / (max_val - min_val);
      }
    }
  }
  
  // Set N k-mer vs N k-mer to 0
  scaled_hamming_matrix(n_kmers - 1, n_kmers - 1) = 0.0;
  
  rownames(scaled_hamming_matrix) = kmer_names;
  colnames(scaled_hamming_matrix) = kmer_names;
  
  return scaled_hamming_matrix;
}