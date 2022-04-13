#!/usr/bin/env Rscript

### Usage: Rscript --vanilla hw1.R <input file> <score file>
### Example: Rscript --vanilla hw1.R input.txt blosum62.txt
### Note: Smith-Waterman Algorithm

### This is one way to read in arguments in R
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("At least two arguments must be supplied (inputFile, scoreFile).n", call.=FALSE) } else if (length(args)>=2) {
    # default gap penalties
    args[3] = -2
    args[4] = -1 }

## Specifying author and email
p <- c(person("Alev", "Baysoy", role = "aut", email = "alev.baysoy@yale.edu"))


# -----------------------------------------------------------------
# runSW
# -----------------------------------------------------------------
#' Title
#'
#' @param inputFile
#' @param scoreFile
#' @param openGap
#' @param extGap
#' @export
runSW = function(inputFile, scoreFile, openGap = -2, extGap = -1) {

  # convert string to integer
  openGap =  strtoi(openGap)
  extGap = strtoi(extGap)

  # Sequences
  cat("\n", file = inputFile, append = TRUE)
  sequences = read.table(file = inputFile, header = FALSE)

  sequence1 = trimws(sequences[1,]) # as string
  sequence2 = trimws(sequences[2,]) # as string

  seq1 = unlist(strsplit(sequence1, "")) # as string vector
  seq2 = unlist(strsplit(sequence2, "")) # as string vector

  # Substitution Matrix
  score_matrix = read.table(file = scoreFile, header = TRUE, sep = "")

  # Construct Matrices
  matrices = construct_matrices(seq1, seq2, score_matrix, openGap, extGap)
  M = matrices$scoring_matrix
  A = matrices$traceback_matrix

  # Traceback
  result = traceback(seq1, seq2, M, A)

  # Write Results
  write2file(filepath = "output.txt",
             sequences = c(sequence1, sequence2),
             scoring_matrix = M,
             sequences_aligned = result)

}

# -----------------------------------------------------------------
# Gap Penalty
# -----------------------------------------------------------------
#' Title
#'
#' @param k
#' @param openGap
#' @param extGap
#' @export
W = function(k, openGap, extGap) {
  return ( (-1) * (openGap + extGap * (k - 1)) )
}

# -----------------------------------------------------------------
# Construct Matrices
# -----------------------------------------------------------------
#' Title
#'
#' @param seq1
#' @param seq2
#' @param score_matrix
#' @param openGap
#' @param extGap
#' @export
#'

construct_matrices = function(seq1, seq2, score_matrix, openGap, extGap) {

  # Dimensions
  col_count = length(seq1) + 1
  row_count = length(seq2) + 1

  # Substitution Matrix (Blosum)
  S = score_matrix

  # Scoring Matrix
  M = matrix(0, nrow = row_count, ncol = col_count)
  M[1,] = 0
  M[,1] = 0
  colnames(M) = c(" ", seq1)
  rownames(M) = c(" ", seq2)

  # Traceback Matrix
  A = matrix(0, nrow = row_count, ncol = col_count)
  A[1,] = "S" # STOP
  A[,1] = "S" # STOP
  colnames(A) = c(" ", seq1)
  rownames(A) = c(" ", seq2)

  # Fill
  for (i in 2:row_count) {
    for (j in 2:col_count){

      parts = integer(4)
      colvec = integer(i)
      rowvec = integer(j)
      for (k in 1:i-1) {colvec[k] = M[i-k,j] - W(k, openGap, extGap)}
      for (k in 1:j-1) {rowvec[k] = M[i,j-k] - W(k, openGap, extGap)}

      parts = c(0,
                M[i-1,j-1] + S[seq2[i-1],seq1[j-1]],
                max(colvec),
                max(rowvec))

      M[i,j] = max(parts)

      if (M[i,j] == parts[1]) {
        A[i,j] = "S" # STOP
      } else if (M[i,j] == parts[2]) {
        A[i,j] = "D" # DIAGONAL
      } else if (M[i,j] == parts[3]) {
        A[i,j] = "U" # UP
      } else if (M[i,j] == parts[4]) {
        A[i,j] = "L" # LEFT
      } else {
        A[i,j] = "E" # ERROR
      }

    }
  }


  return (list("scoring_matrix" = M, "traceback_matrix" = A))

}

# -----------------------------------------------------------------
# Traceback
# -----------------------------------------------------------------
#' Title
#'
#' @param seq1
#' @param seq2
#' @param M
#' @param A
#' @export
traceback = function(seq1, seq2, M, A) {

  # locate maximum value
  z = which(M == max(M), arr.ind = TRUE)
  i = z[1]
  j = z[2]

  s1 = ""
  s2 = ""
  bars = ""
  space = ""

  while (A[i, j] != "S") {

    if (A[i,j] == "D") {
      s2 = paste(seq2[i-1], s2, sep = space)
      s1 = paste(seq1[j-1], s1, sep = space)
      if (seq2[i-1] == seq1[j-1]) {
        bars = paste("|", bars, sep = space)
      } else {
        bars = paste(" ", bars, sep = space)
      }
      i = i-1
      j = j-1

    } else if (A[i,j] == "U") {
      s2 = paste(seq2[i-1], s2, sep = space)
      s1 = paste("-", s1, sep = space)
      bars = paste(" ", bars, sep = space)
      i = i-1

    } else if (A[i,j] == "L") {
      s2 = paste("-", s2, sep = space)
      s1 = paste(seq1[j-1], s1, sep = space)
      bars = paste(" ", bars, sep = space)
      j = j-1

    } else if (A[i,j] == "E") {
      stop("Traceback matrix contains an erroneous value.")
    }

  }

  return (c(s1, s2, bars))

}

# -----------------------------------------------------------------
# Write Results to File
# -----------------------------------------------------------------
#' Title
#'
#' @param filepath
#' @param sequences
#' @param scoring_matrix
#' @param sequences_aligned
#' @export
#'
write2file = function(filepath, sequences, scoring_matrix, sequences_aligned) {

  # Start writing to an output file
  sink(filepath)

  cat("-----------", sep = "\n")
  cat("|Sequences|", sep = "\n")
  cat("-----------", sep = "\n")
  cat("sequence1", sep = "\n")
  cat(sequences[1], sep = "\n")
  cat("sequence2", sep = "\n")
  cat(sequences[2], sep = "\n")

  cat("--------------", sep = "\n")
  cat("|Score Matrix|", sep = "\n")
  cat("--------------", sep = "\n")
  data = scoring_matrix
  cat(" ", colnames(data), sep = "\t")
  cat("\n")
  cat(write.table(data, file = "", append = TRUE, sep = "\t", row.names=T, col.names=F, quote=FALSE))

  cat("----------------------", sep = "\n")
  cat("|Best Local Alignment|", sep = "\n")
  cat("----------------------", sep = "\n")
  cat("Alignment Score:")
  cat(max(scoring_matrix), sep = "\n")
  cat("Alignment Results: ", sep = "\n")
  formatted = format(sequences, sequences_aligned)
  cat(formatted$sequence1, sep = "\n")
  cat(formatted$bars, sep = "\n")
  cat(formatted$sequence2, sep = "\n")

  # Stop writing to the file
  sink()

}

# -----------------------------------------------------------------
# Format
# -----------------------------------------------------------------
#' Title
#'
#' @param sequences
#' @param sequences_aligned
#'
#' @export

format = function(sequences, sequences_aligned) {

  # reconstruct sequence1 with insertions added
  s1 = gsub("-", "", sequences_aligned[1])
  substring = paste("(", sequences_aligned[1], ")", sep = "")
  s1 = gsub(s1, substring, sequences[1])

  # reconstruct sequence2 with insertions added
  s2 = gsub("-", "", sequences_aligned[2])
  substring = paste("(", sequences_aligned[2], ")", sep = "")
  s2 = gsub(s2, substring, sequences[2])

  # original bar sequence from traceback
  bars = sequences_aligned[3]

  # ---------------------------------------
  # add padding to align sequences at "("
  # ---------------------------------------

  L1 = which(strsplit(s1, "")[[1]] == "(")   # position of "(" in s1
  L2 = which(strsplit(s2, "")[[1]] == "(")   # position of "(" in s2
  padding_length = abs(L1 - L2)

  if (L1 >= L2) {
    if (L2 > 0) {for (k in 1:L2) {bars = paste(" ", bars, sep = "")}}
    if (padding_length > 0) {
      for (k in 1:padding_length) {
        s2 = paste(" ", s2, sep = "")
        bars = paste(" ", bars, sep = "")
      }
    }

  } else {
    if (L1 > 0) {for (k in 1:L1) {bars = paste(" ", bars, sep = "")}}
    if (padding_length > 0) {
      for (k in 1:padding_length) {
        s1 = paste(" ", s1, sep = "")
        bars = paste(" ", bars, sep = "")
      }
    }
  }

  return (list("sequence1" = s1, "sequence2" = s2, "bars" = bars))

}


# -----------------------------------------------------------------
# Main
# -----------------------------------------------------------------
runSW(inputFile=args[1], scoreFile=args[2], openGap=args[3], extGap=args[4])


