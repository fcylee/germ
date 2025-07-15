#include <Rcpp.h>
using namespace Rcpp;

//' Chops sequence into k-mers
//'
//' @param input_seq a character string
//' @param k_len an integer specifying the length of the k-mer
//'
//' @return character vector of k-mers
//' @export
// [[Rcpp::export]]
CharacterVector kmer_chopper(std::string input_seq, int k_len)
{
  int num_kmers = input_seq.length() - (k_len - 1);
  CharacterVector output_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    output_kmers[i] = input_seq.substr (i, k_len);
  }

  return(output_kmers);
}

//' Calculates k-mer multivalencies
//'
//' @param input_seq sequence string
//' @param k_len an integer specifying the length of the k-mer
//' @param window_size integer specifying window_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance vector
//'
//' @return a vector of k-mer multivalencies
//' @export
// [[Rcpp::export]]
NumericVector calculate_kmer_multivalencies(std::string input_seq, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances)
{
  // First, split the input string up into k-mers, then convert those k-mers into their indices on the Hamming distance matrix.

  int num_kmers = input_seq.length() - (k_len - 1);
  CharacterVector input_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    input_kmers[i] = input_seq.substr (i, k_len);
  }

  // CharacterVector input_kmers = kmer_chopper2(input_seq, k_len);
  CharacterVector hd_kmers = rownames(hamming_distances);
  NumericVector match_kmers(input_kmers.size());

  match_kmers = match(input_kmers, hd_kmers);

  // If there are unexpected characters in the input sequence, the match function returns a negative value.
  // We can set these values to the index of the N kmer, so that they will be counted as 0 when we score the multivalency.
  int hd_dim = pow(4, k_len) + 1; // Dimensions of the Hamming distance matrix (including N kmer)

  for (int i = 0; i < num_kmers; ++i) {
  	if (match_kmers[i] < 0) {
  	match_kmers[i] = hd_dim;
  	}
  }


  // Remember that R indexes from 1, not 0, so the input_kmers are 1 larger than their index in the Hamming distance matrix.

  // Initialise an output vector. The output vector should be the size of the input k-mers before N-padding.
  NumericVector output_vector(input_kmers.size());

  // Initialise all the variables
  int central_kmer = round((window_size - 0.1) / 2);
  int central_kmer_index;
  double sum_score;

  NumericVector local_kmers(window_size);
  double ham_value;
  double pos_value;
  double score;

  // Pad the input kmer indices with the index that corresponds to the N-kmer (to permit calculation of edges).
  // int hd_dim = pow(4, k_len) + 1;

  NumericVector padded_kmers(match_kmers.size() + (2 * central_kmer));

  // Shockingly, this seems to be the fastest way to concatenate two vectors.
  for(int i = 0; i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  for(int i = central_kmer; i < central_kmer + match_kmers.size(); ++i) {
    padded_kmers[i] = match_kmers[i - central_kmer];
  }

  for(int i = central_kmer + match_kmers.size(); i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  // Here is where the actual work happens.
  for(int i = 0; i < output_vector.size(); ++i) {
    // Take a window of indices.
    local_kmers = padded_kmers[Rcpp::Range(i, i + window_size - 1)];
    // Get the central index.
    central_kmer_index = local_kmers[central_kmer] - 1;

    sum_score = 0;

    for(int j = 0; j < window_size; ++j) {
      // Get the hamming distances between every k-mer and the central k-mer. Multiply by the positional distance. Add the result to the scoring vector.
      ham_value = hamming_distances(local_kmers[j] - 1, central_kmer_index);
      pos_value = positional_distances[j];
      score =  ham_value * pos_value;
      sum_score += score;
    }

    output_vector[i] = sum_score;
  }

  return(output_vector);
}

//' Get list of k-mer multivalencies
//'
//' @param ins list of ??
//' @param k_len an integer specifying the length of the k-mer
//' @param window_size integer specifying window_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance vector
//'
//' @return a list of k-mer multivalencies
//' @export
// [[Rcpp::export]]
List list_kmer_multivalencies(List ins, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances)
{
  // Dedicated replacement for lapply. This is faster when using a larger number of input sequences.

  List output_list(ins.size());

  for(int i = 0; i < ins.size(); ++i){

    output_list[i] = calculate_kmer_multivalencies(ins[i], k_len, window_size, hamming_distances, positional_distances);

  }

  return(output_list);

}

//' Calculates sliding mean and pads ends with 0
//'
//' @param iv numeric vector
//' @param ws sliding window size
//'
//' @return numeric vector of sliding means padded with 0
//' @export
// [[Rcpp::export]]
NumericVector calculate_padded_sliding_mean(NumericVector iv, int ws)
{
  // If smoothing window is less than length of numeric vector return all NA
  if(iv.size() < ws) {

    NumericVector pad(iv.size());
    for(int i = 0; i < iv.size(); ++i) {
      pad[i] = NA_REAL;
    }
    return pad;

  } else {
  
    // Simple sliding window function to compute means of vector iv in windows of size ws. Returns vector that is ws - 1 shorter than the input vector.
    int n = iv.size();
    NumericVector out(n - (ws - 1));

    NumericVector wv = iv[seq(0, ws - 1)];
    double sv = sum(wv);
    out[0] = sv/ws;

    for(int i = 1; i < n - (ws - 1); ++i) {
      sv += iv[i + ws - 1] - iv[i-1];
      out[i] = sv/ws;
    }

    NumericVector pad(iv.size());
    for(int i = (ws - 1)/2; i < iv.size() - (ws - 1)/2; ++i) {
      pad[i] = out[i - (ws - 1)/2];
    }
    return pad;

  }

}

//' Calculates k-mer multivalencies with tidy output
//'
//' @param input_seq sequence string
//' @param input_seq_name sequence string name (e.g. transcript id)
//' @param k_len an integer specifying the length of the k-mer
//' @param window_size integer specifying window_size
//' @param smoothing size integer specifying smoothingwindow_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance vector
//'
//' @return a data frame of k-mer multivalencies
//' @export
// [[Rcpp::export]]
DataFrame calculate_kmer_multivalencies_df(std::string input_seq, std::string input_seq_name, int k_len, int window_size, int smoothing_size, NumericMatrix hamming_distances, NumericVector positional_distances)
{
  // First, split the input string up into k-mers, then convert those k-mers into their indices on the Hamming distance matrix.

  int num_kmers = input_seq.length() - (k_len - 1);
  CharacterVector input_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    input_kmers[i] = input_seq.substr (i, k_len);
  }

  // CharacterVector input_kmers = kmer_chopper2(input_seq, k_len);
  CharacterVector hd_kmers = rownames(hamming_distances);
  NumericVector match_kmers(input_kmers.size());

  match_kmers = match(input_kmers, hd_kmers);

  // If there are unexpected characters in the input sequence, the match function returns a negative value.
  // We can set these values to the index of the N kmer, so that they will be counted as 0 when we score the multivalency.
  int hd_dim = pow(4, k_len) + 1; // Dimensions of the Hamming distance matrix (including N kmer)

  for (int i = 0; i < num_kmers; ++i) {
  	if (match_kmers[i] < 0) {
  	match_kmers[i] = hd_dim;
  	}
  }


  // Remember that R indexes from 1, not 0, so the input_kmers are 1 larger than their index in the Hamming distance matrix.

  // Initialise an output vector. The output vector should be the size of the input k-mers before N-padding.
  NumericVector output_vector(input_kmers.size());

  // Initialise all the variables
  int central_kmer = round((window_size - 0.1) / 2);
  int central_kmer_index;
  double sum_score;

  NumericVector local_kmers(window_size);
  double ham_value;
  double pos_value;
  double score;

  // Pad the input kmer indices with the index that corresponds to the N-kmer (to permit calculation of edges).
  // int hd_dim = pow(4, k_len) + 1;

  NumericVector padded_kmers(match_kmers.size() + (2 * central_kmer));

  // Shockingly, this seems to be the fastest way to concatenate two vectors.
  for(int i = 0; i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  for(int i = central_kmer; i < central_kmer + match_kmers.size(); ++i) {
    padded_kmers[i] = match_kmers[i - central_kmer];
  }

  for(int i = central_kmer + match_kmers.size(); i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  // Here is where the actual work happens.
  for(int i = 0; i < output_vector.size(); ++i) {
    // Take a window of indices.
    local_kmers = padded_kmers[Rcpp::Range(i, i + window_size - 1)];
    // Get the central index.
    central_kmer_index = local_kmers[central_kmer] - 1;

    sum_score = 0;

    for(int j = 0; j < window_size; ++j) {
      // Get the hamming distances between every k-mer and the central k-mer. Multiply by the positional distance. Add the result to the scoring vector.
      ham_value = hamming_distances(local_kmers[j] - 1, central_kmer_index);
      pos_value = positional_distances[j];
      score =  ham_value * pos_value;
      sum_score += score;
    }

    output_vector[i] = sum_score;
  }

  // Now tidy into a data frame output
  CharacterVector seqnames(input_kmers.size());
  for(int i = 0; i < input_kmers.size(); ++i) {
    seqnames[i] = input_seq_name;
  }

  // Calculate sliding mean of scores
  NumericVector smoothed_vector = input_kmers.size();
  smoothed_vector = calculate_padded_sliding_mean(output_vector, smoothing_size);

  DataFrame df = DataFrame::create(Named("sequence_name") = seqnames,
                                   Named("kmer") = input_kmers,
                                   Named("kmer_multivalency") = output_vector,
                                   Named("smoothed_kmer_multivalency") = smoothed_vector);

  return(df);
}

//' Calculates k-mer multivalencies with tidy output for complete 2d matrix and a positional distance matrix
//'
//' @param input_seq sequence string
//' @param input_seq_name sequence string name (e.g. transcript id)
//' @param k_len an integer specifying the length of the k-mer
//' @param smoothing size integer specifying smoothingwindow_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance matrix
//'
//' @return a list with matrix of pairwise k-mer weighted similarities and data frame of k-mer multivalencies summed per position
//' @export
// [[Rcpp::export]]
List calculate_kmer_pairwise_multivalencies(std::string input_seq, std::string input_seq_name, int k_len, int smoothing_size, NumericMatrix hamming_distances, NumericMatrix positional_distances)
{
  // Calculate centering offset
  int center_offset = (k_len - 1) / 2;
  int seq_length = input_seq.length();
  
  // Extract k-mers (standard extraction for mapping)
  int num_kmers = seq_length - (k_len - 1);
  CharacterVector input_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    input_kmers[i] = input_seq.substr(i, k_len);
  }

  // Map k-mers to indices
  CharacterVector hd_kmers = rownames(hamming_distances);
  NumericVector match_kmers(input_kmers.size());
  match_kmers = match(input_kmers, hd_kmers);

  // Handle unexpected characters
  int hd_dim = pow(4, k_len) + 1;
  for (int i = 0; i < num_kmers; ++i) {
    if (match_kmers[i] < 0) {
      match_kmers[i] = hd_dim;
    }
  }

  // Create full sequence length pairwise matrix
  NumericMatrix output_matrix(seq_length, seq_length);

  for(int i = 0; i < seq_length; ++i) {
    for(int j = 0; j < seq_length; ++j) {
      
      // Determine k-mer indices for positions i and j
      int kmer_index_i, kmer_index_j;
      
      // Check if position i can have a centered k-mer
      if (i >= center_offset && i < center_offset + num_kmers) {
        int kmer_pos_i = i - center_offset;  // Which k-mer is centered at position i
        kmer_index_i = match_kmers[kmer_pos_i] - 1;
      } else {
        kmer_index_i = hd_dim - 1;  // N k-mer (should be zeros)
      }
      
      // Same for position j
      if (j >= center_offset && j < center_offset + num_kmers) {
        int kmer_pos_j = j - center_offset;
        kmer_index_j = match_kmers[kmer_pos_j] - 1;
      } else {
        kmer_index_j = hd_dim - 1;  // N k-mer (should be zeros)
      }
      
      // Calculate similarity
      double ham_value = hamming_distances(kmer_index_i, kmer_index_j);
      double pos_weight = positional_distances(i, j);
      
      output_matrix(i, j) = ham_value * pos_weight;
    }
  }

  // Calculate summed scores per position
  NumericVector position_scores(seq_length);
  for(int i = 0; i < seq_length; ++i) {
    double sum_score = 0;
    for(int j = 0; j < seq_length; ++j) {
      sum_score += output_matrix(i, j);
    }
    position_scores[i] = sum_score;
  }

  // Calculate sliding mean of scores (same as original)
  NumericVector smoothed_scores(seq_length);
  smoothed_scores = calculate_padded_sliding_mean(position_scores, smoothing_size);

  // Create k-mer column for each position (like original)
  CharacterVector seqnames_pos(seq_length);
  CharacterVector position_kmers(seq_length);
  std::string n_kmer(k_len, 'N');  // N k-mer for edge positions
  
  for(int i = 0; i < seq_length; ++i) {
    seqnames_pos[i] = input_seq_name;
    
    // Get k-mer centered at position i
    if (i >= center_offset && i < center_offset + num_kmers) {
      int kmer_pos = i - center_offset;
      position_kmers[i] = input_kmers[kmer_pos];
    } else {
      position_kmers[i] = n_kmer;  // N k-mer for edges
    }
  }

  // Return position summary with k-mers (like original code)
  DataFrame df = DataFrame::create(Named("sequence_name") = seqnames_pos,
                                   Named("kmer") = position_kmers,
                                   Named("position_multivalency") = position_scores,
                                   Named("smoothed_position_multivalency") = smoothed_scores);

  List result = List::create(Named("matrix") = output_matrix,
                            Named("position_summary") = df);
return(result);
}
