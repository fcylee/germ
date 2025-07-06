#' Plots smoothed k-mer multivalency for a transcript
#'
#' @param kmer_multivalency.dt k-mer multivalency dataframe/data.table
#' @param k_len k-mer length
#' @param seq Character vector of sequences from FASTA
#' @param seq_name Sequence name
#' @param outdir Output plot directory
#' @param interactive If TRUE, will also plot to Rstudio window
#'
#' @return PDF of smoothed k-mer multivalency plot
#' @export
#' @import ggplot2
plot_kmer_multivalency_pdm <- function(kmer_multivalency.dt, k_len, seq, seq_name, outdir, interactive = FALSE, annotate_max = FALSE) {
  tx.dt <- kmer_multivalency.dt[kmer_multivalency.dt$sequence_name == seq_name, ]
  
  # CHANGED: Removed offset since new approach returns one row per sequence position, not per k-mer
  tx.dt$coord <- 1:nrow(tx.dt)  # Was: 1:nrow(tx.dt) + floor(k_len/2)
  
  # CHANGED: Simplified sanity check since we now have seq_length rows instead of num_kmers rows
  stopifnot(max(tx.dt$coord) == nchar(seq[names(seq) == seq_name])) # Was: (max(tx.dt$coord) + floor(k_len/2)) == nchar(seq[names(seq) == seq_name])
  
  p <- ggplot() +
    # CHANGED: Updated column name to match new DataFrame output
    geom_line(data = tx.dt, aes(x = coord, y = smoothed_position_multivalency)) +  # Was: smoothed_kmer_multivalency
    labs(x = seq_name,
         y = paste0(k_len, "-mer multivalency")) +
    coord_cartesian(xlim = c(1, nchar(seq[names(seq) == seq_name]))) +
    theme_minimal()
  
  if(annotate_max) {
    p <- p +
      # CHANGED: Updated column name to match new DataFrame output
      geom_vline(xintercept = tx.dt[which.max(tx.dt$smoothed_position_multivalency), ]$coord,  # Was: smoothed_kmer_multivalency
                 linetype = "dashed", colour = "grey50") +
      # CHANGED: Updated column name to match new DataFrame output
      geom_text(data = tx.dt[which.max(tx.dt$smoothed_position_multivalency), ],  # Was: smoothed_kmer_multivalency
                aes(x = coord, y = 0, label = kmer),
                angle = 90, hjust = 0, vjust = 1.5)
  }
  
  ggsave(plot = p,
         filename = paste0(outdir, "/", seq_name, ".pdf"),
         width = 200,
         height = 75,
         units = "mm")
  if(interactive) print(p)
}