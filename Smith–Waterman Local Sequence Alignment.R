############################################################
# Smith–Waterman Local Sequence Alignment
# Author: Jordan Scollick-Karnitz
# Language: R
#
# Description:
# Implements the Smith–Waterman algorithm for local
# sequence alignment using dynamic programming.
# Supports multiple optimal local alignments.
############################################################

local_alignment <- function(seq1 = "ATCG",
                            seq2 = "ATG",
                            match = 3,
                            mismatch = -3,
                            gap = -2) {
  
  # Sequence lengths
  n <- nchar(seq1)
  m <- nchar(seq2)
  
  # Split sequences into characters
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  # Initialize score matrix
  F <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  # Traceback matrix
  traceback <- vector("list", (n + 1) * (m + 1))
  dim(traceback) <- c(n + 1, m + 1)
  
  # Track maximum score(s)
  max_score <- 0
  max_positions <- list()
  
  # Fill scoring and traceback matrices
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      
      score_diag <- F[i - 1, j - 1] +
        ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
      
      score_up   <- F[i - 1, j] + gap
      score_left <- F[i, j - 1] + gap
      
      best_score <- max(0, score_diag, score_up, score_left)
      F[i, j] <- best_score
      
      directions <- character(0)
      if (best_score > 0) {
        if (score_diag == best_score) directions <- c(directions, "diag")
        if (score_up   == best_score) directions <- c(directions, "up")
        if (score_left == best_score) directions <- c(directions, "left")
      }
      
      traceback[[i, j]] <- directions
      
      # Track positions with maximum score
      if (best_score > max_score) {
        max_score <- best_score
        max_positions <- list(c(i, j))
      } else if (best_score == max_score && best_score > 0) {
        max_positions <- append(max_positions, list(c(i, j)))
      }
    }
  }
  
  alignments <- list()
  
  # Recursive traceback for local alignments
  trace <- function(i, j, aln1 = "", aln2 = "") {
    
    # Stop traceback when score reaches zero
    if (F[i, j] == 0) {
      alignments <<- append(alignments, list(c(aln1, aln2)))
      return()
    }
    
    for (d in traceback[[i, j]]) {
      if (d == "diag") {
        trace(i - 1, j - 1,
              paste0(s1[i - 1], aln1),
              paste0(s2[j - 1], aln2))
      } else if (d == "up") {
        trace(i - 1, j,
              paste0(s1[i - 1], aln1),
              paste0("-", aln2))
      } else if (d == "left") {
        trace(i, j - 1,
              paste0("-", aln1),
              paste0(s2[j - 1], aln2))
      }
    }
  }
  
  # Trace back from all maximum-scoring cells
  for (pos in max_positions) {
    trace(pos[1], pos[2])
  }
  
  return(list(
    score = max_score,
    alignments = alignments,
    score_matrix = F
  ))
}

############################################################
# Example usage
############################################################

cat("LOCAL SEQUENCE ALIGNMENT (Smith–Waterman)\n\n")

result <- local_alignment("ATCG", "ATG")

cat("Sequence 1: ATCG\n")
cat("Sequence 2: ATG\n\n")
cat("Best Local Alignment Score:", result$score, "\n\n")
cat("Number of optimal local alignments:", length(result$alignments), "\n\n")

for (i in seq_along(result$alignments)) {
  aln <- result$alignments[[i]]
  cat("Alignment", i, ":\n")
  cat("Seq1:", aln[1], "\n")
  cat("Seq2:", aln[2], "\n")
  
  matches <- ifelse(
    strsplit(aln[1], "")[[1]] ==
      strsplit(aln[2], "")[[1]],
    "|", " "
  )
  cat("      ", paste(matches, collapse = ""), "\n\n")
}

cat("Score Matrix:\n")
print(result$score_matrix)
