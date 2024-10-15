CC = read.csv("peaklists/CC_empty.csv")
NCACX = read.csv("peaklists/NCACX_empty.csv")

print(paste("Total number of NCACX peaks:", nrow(NCACX)))

is_adj = function(v1, v2, tolerance=1) {
  # return true if the vectors are neighboring, false otherwise
  condition_1 = all(abs(v1 - v2) <= tolerance)
  condition_2 = all(abs(v1[length(v1):1] - v2) <= tolerance)

  if (condition_1 | condition_2) {
    # if the locations match in any order of x and y, return TRUE
    return (TRUE)
  }
  return (FALSE) 
}

map_to_CC = function(v1, tolerance=1) {
  for (row_2 in 1:nrow(CC)) {
    v2 = CC[row_2, c("Position.F1", "Position.F2")]
    if (is_adj(v1, v2, tolerance)) {
      print(paste(TRUE, v1, v2))
      return (TRUE)
    }
  }
  return (FALSE)
}

# is_adj(v1=c(1, 2), v2=c(1, 1))
# is_adj(v1=c(1, 2), v2=c(2, 8))

map_NCACX_to_CC = function(tolerance=1) {
  # return a vector of booleans indicating if each NCACX peak exist in CC

  out = c()
  for (row_1 in 1:nrow(NCACX)) {
    print(paste("row_1=", row_1))
    v1 = NCACX[row_1, c("Position.F1", "Position.F2")]
    is_in_CC = map_to_CC(v1, tolerance)
    print(paste("is_in_CC=", is_in_CC))
    out = append(out, map_to_CC(v1, tolerance))
  }
  return (out)
}


filter_vec = map_NCACX_to_CC(tolerance=1)
subset = subset(NCACX, filter_vec)
print(paste("Number of peaks remaining after filtering:", nrow(subset)))
      
write.csv(subset, "NCACX_processed.csv")
print("Filtered data saved to csv.")
