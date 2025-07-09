library(reticulate)
library(Biostrings)
library(ggplot2)     # for ggsave
np <- import("numpy")

# takes a vector of logits (real #'s) and transforms into prob distr
softmax <- function(x) {
  # subtracting each element from its max (for numerical stability)
  exp_x <- exp(x - max(x))
  exp_x / sum(exp_x) # normalizing vector to sum to 1
}

get_logits <- function(filename) {
  base_np <- np$load(filename) # loading the function and putting into base_np

  # convert numpy file to Rx
  mat <- py_to_r(base_np)[1, , ] # python to r array (takes first slice)

  # turn into probability matrix using softmax (N samples x M tokens)
  # applies 3rd arg to 1st arg (1 = row wise) & t is transpose
  all <- t(apply(mat, 1, softmax))
  # focus only on nucleotides 
  nuc_prob <- data.frame(
    A = all[, utf8ToInt("A") + 1], # utf8ToInt ASCII + of "A" and then adding 1
    C = all[, utf8ToInt("C") + 1],
    G = all[, utf8ToInt("G") + 1],
    T = all[, utf8ToInt("T") + 1]
  )

  # shifting EVO2 data: EVO2 data starts at pos 2 (FASTA data starts at pos 1)
  empty_row <- data.frame(A = 0.25, C = 0.25, G = 0.25, T = 0.25) # at pos 1
  nuc_prob <- rbind(empty_row, nuc_prob)
  nuc_prob[-nrow(nuc_prob), ]
}

#' Save a ggplot object to PDF in the specified folder (creates it if absent)
save_as_pdf <- function(p, filename, out_dir = "My_figures", width = 8, height = 6) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggsave(file.path(out_dir, filename), plot = p, device = "pdf", width = width, height = height)
}