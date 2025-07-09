#' Plots k-mer multivalency matrix as heatmap with mean/max split by diagonal
#'
#' @param kmer_multivalency_all List output from calculate_kmer_pairwise_multivalencies
#' @param seq_name Sequence name to plot
#' @param k_len K-mer length
#' @param outdir Output plot directory
#' @param interactive If TRUE, will also plot to Rstudio window
#' @param max_size Maximum matrix size to plot before binning
#'
#' @return PNG of k-mer multivalency heatmap (if binned - mean upper, max lower)
#' @export
#' @import ggplot2
#' @import viridis
plot_kmer_multivalency_pdm_heatmap <- function(kmer_multivalency_all, seq_name, k_len, outdir, 
                                               interactive = FALSE, max_size = 2000, percentile_clip = c(0.10, 0.99)) {
  
  # Extract matrix from result
  multivalency_matrix <- kmer_multivalency_all$matrix
  seq_length <- nrow(multivalency_matrix)
  original_size <- seq_length  # Store original size
  binning_done <- FALSE
  bin_size <- 1
  
  # Bin if matrix is too large
  if(seq_length > max_size) {
    binning_done <- TRUE
    # Calculate bin size
    bin_size <- ceiling(seq_length / max_size)
    
    # Create binned matrices for both methods
    new_size <- floor(seq_length / bin_size)
    binned_matrix_mean <- matrix(0, nrow = new_size, ncol = new_size)
    binned_matrix_max <- matrix(0, nrow = new_size, ncol = new_size)
    
    for(i in 1:new_size) {
      for(j in 1:new_size) {
        # Define bin boundaries
        row_start <- (i-1) * bin_size + 1
        row_end <- min(i * bin_size, seq_length)
        col_start <- (j-1) * bin_size + 1
        col_end <- min(j * bin_size, seq_length)
        
        # Apply both binning methods
        bin_values <- multivalency_matrix[row_start:row_end, col_start:col_end]
        binned_matrix_mean[i, j] <- mean(bin_values)
        binned_matrix_max[i, j] <- max(bin_values)
      }
    }
    
    message(paste("Matrix binned from", seq_length, "x", seq_length, 
                  "to", new_size, "x", new_size, "using", bin_size, "x", bin_size, "bins"))
    
    # Update sequence length
    seq_length <- new_size
  } else {
    # No binning needed - use original matrix for both
    binned_matrix_mean <- multivalency_matrix
    binned_matrix_max <- multivalency_matrix
  }
  
  # Create combined matrix: upper triangle = mean, lower triangle = max
  combined_matrix <- matrix(0, nrow = seq_length, ncol = seq_length)
  
  for(i in 1:seq_length) {
    for(j in 1:seq_length) {
      if(i <= j) {
        combined_matrix[i, j] <- binned_matrix_mean[i, j]  # Upper triangle (including diagonal)
      } else {
        combined_matrix[i, j] <- binned_matrix_max[i, j]   # Lower triangle
      }
    }
  }
  
  # Convert to long format for ggplot
  heatmap_data <- data.frame(
    position_i = rep(1:seq_length, seq_length),
    position_j = rep(1:seq_length, each = seq_length),
    multivalency = as.vector(combined_matrix)
  )
  
  # Create subtitle with binning information
  if(binning_done) {
    subtitle_text <- paste("Upper triangle: Mean | Lower triangle: Max | Binned:", 
                          original_size, "→", seq_length, "(", bin_size, "×", bin_size, "bins )")
  } else {
    subtitle_text <- "Upper triangle: Mean | Lower triangle: Max | No binning"
  }
  
  # Create the heatmap
  p <- ggplot(heatmap_data, aes(x = position_i, y = position_j, fill = multivalency)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = paste0("Distance-Weighted ", k_len, "-mer Similarity: \n", seq_name),
         subtitle = subtitle_text,
         x = "Sequence/Binned Position",
         y = "Sequence/Binned Position",
         fill = "Similarity\nScore") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed()
  
  # Calculate appropriate plot dimensions
  plot_width <- min(12, max(6, seq_length * 0.002))
  plot_height <- plot_width  # Square
  dpi <- min(300, max(150, 1500 / seq_length))
  
  # Save as PNG
  ggsave(plot = p,
         filename = paste0(outdir, "/", seq_name, "_heatmap.png"),
         width = plot_width,
         height = plot_height,
         units = "in",
         dpi = dpi)
         
  # Create clipped version with percentile-based color scale
  percentiles <- quantile(heatmap_data$multivalency, probs = percentile_clip, na.rm = TRUE)
  clipped_min <- percentiles[1]
  clipped_max <- percentiles[2]

  # Create the clipped heatmap
  p_clipped <- ggplot(heatmap_data, aes(x = position_i, y = position_j, fill = multivalency)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(clipped_min, clipped_max)) +  # Use percentile limits
    labs(title = paste0("Distance-Weighted ", k_len, "-mer Similarity: \n", seq_name, " (Clipped)"),
         subtitle = paste0(subtitle_text, "\nColor scale: ", min(percentile_clip)*100, "th-", max(percentile_clip)*100, "th percentile"),
         x = "Sequence/Binned Position",
         y = "Sequence/Binned Position",
         fill = "Similarity\nScore") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    coord_fixed()

  # Save clipped version
  ggsave(plot = p_clipped,
         filename = paste0(outdir, "/", seq_name, "_heatmap_clipped.png"),
         width = plot_width,
         height = plot_height,
         units = "in",
         dpi = dpi)

  message(paste("Clipped color scale range:", round(clipped_min, 3), "to", round(clipped_max, 3)))
  
  if(interactive) print(p)
  
  return(p)
}