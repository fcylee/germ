#' Creates positional distance matrix based on window size
#'
#' @param window_size an integer specifying the window size
#' @param k_len an integer specifying the length of the k-mer
#' @return scaled distance integer matrix
#' @export
#'
create_positional_distance_matrix <- function(distance_matrix, window_size, k_len) {
  seq_length <- nrow(distance_matrix)
  
  # Apply triangular scaling based on physical window
  weight_matrix <- (window_size - distance_matrix) / window_size
  weight_matrix[weight_matrix < 0] <- 0  # Clip distances > window_size
  
  # Zero out overlapping k-mers based on sequence position difference
  seq_positions <- 1:seq_length
  position_diff_matrix <- abs(outer(seq_positions, seq_positions, "-"))
  weight_matrix[position_diff_matrix < k_len] <- 0
  
  return(weight_matrix)
}