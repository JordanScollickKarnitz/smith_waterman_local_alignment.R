# Smith–Waterman Local Sequence Alignment (R)

This repository contains an implementation of the **Smith–Waterman algorithm**
for local sequence alignment written in **R**.

---

## 🧬 Features
- Local sequence alignment
- Customizable scoring scheme
- Supports multiple optimal alignments
- Returns best alignment score and scoring matrix
- Designed for instructional use

---

## 🚀 Usage

```r
source("smith_waterman_local_alignment.R")

result <- local_alignment("ATCG", "ATG")
result$score
result$alignments
result$score_matrix

## ⏱️ Runtime & Space Complexity

Let **n** be the length of sequence 1 and **m** be the length of sequence 2.

### Time Complexity
- The dynamic programming matrix is filled by evaluating three scores
  (diagonal, up, left) for each cell.
- This results in a time complexity of **O(n × m)**.

- Traceback begins from the highest-scoring cell(s) and stops when a cell with
  score zero is reached.
  - In the worst case, multiple optimal local alignments may exist.
  - However, traceback is typically much faster than matrix construction.

### Space Complexity
- The scoring matrix requires **O(n × m)** space.
- The traceback matrix also requires **O(n × m)** space.
- Total space complexity is **O(n × m)**.

### Practical Notes
- Local alignment identifies the most similar subsequences rather than forcing
  alignment across full sequence lengths.
- This makes Smith–Waterman particularly useful for detecting conserved regions
  or motifs within otherwise dissimilar sequences.
