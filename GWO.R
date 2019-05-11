GWO <- function(FUN, optimType = "MIN", numVar, rangeVar, control = list(), numBest = 3) {
  t0 <- Sys.time()

  # Check input variables
  if (!optimType %in% c("MIN", "MAX"))
    error("Parameter 'optimType' must be 'MIN' or 'MAX'")

  if (!is.numeric(numVar))
    error("Parameter 'numVar' must be number")

  if (any(dim(rangeVar) != c(2, numVar)))
    error("Parameter 'rangeVar' must be matrix 2 x numVar")

  if (!is.numeric(rangeVar))
    error("Parameter 'rangeVar' must be matrix of numbers")

  if (!"numPopulation" %in% names(control))
    control$numPopulation <- 40

  if (!"maxIter" %in% names(control))
    control$maxIter <- 500

  if (numBest > control$numPopulation)
    error("Parameter 'numPopulation' must be greater than 'numBest'")

  # Initialization
  desc <- optimType != "MIN"
  i <- 0
  bar <- txtProgressBar(min = 0, max = control$maxIter, style = 3)
  verbose <- rep(0, control$maxIter - 1)

  # Initialization of population
  wolves <- apply(rangeVar, 2, function(x) runif(control$numPopulation, x[1], x[2]))

  # Calculate fitness for each wolf
  solutions <- apply(wolves, 1, function(x) FUN(x))

  # Find best wolves
  wolves <- cbind(wolves, solutions)
  dimnames(wolves) <- NULL
  wolves <- wolves[order(wolves[, numVar + 1], decreasing = desc), ]
  bestWolves <- wolves[1:numBest, ]

  # Compute best wolves in each iteration and update vectors
  while (i < control$maxIter) {
    # Update coefficient
    a <- 2 - (i * 2 / control$maxIter)

    # Update wolves position
    for (j in 1:numVar) {
      # Generate random vectors
      r1 <- matrix(runif(control$numPopulation * numBest, 0, 1), nrow = control$numPopulation, ncol = numBest)
      r2 <- matrix(runif(control$numPopulation * numBest, 0, 1), nrow = control$numPopulation, ncol = numBest)

      # Update coefficients
      A <- 2 * a * r1 - a
      C <- 2 * r2

      # Update wolves position based on best positions
      D <- abs(C * bestWolves[, j] - wolves[, j]) # warning
      X <- bestWolves[, j] - A * D # warning
      wolves[, j] <- apply(X, 1, function(x) sum(x)/length(x))

      # Check boundaries
      wolves[(wolves[, j] < rangeVar[1, j]), j] <- rangeVar[1, j]
      wolves[(wolves[, j] > rangeVar[2, j]), j] <- rangeVar[2, j]
    }

    # Calculate fitness for each wolf
    wolves[, numVar + 1] <- apply(wolves[, 1:numVar], 1, function(x) FUN(x))

    # Find best wolves
    wolves <- wolves[order(wolves[, numVar + 1], decreasing = desc), ]
    bestWolves <- wolves[1:numBest, ]
    verbose[i] <- wolves[1, numVar + 1]

    i <- i + 1
    setTxtProgressBar(bar, i)
  }

  close(bar)
  return(list(
    verbose = verbose,
    result = wolves[1, 1:numVar],
    optimumValue = wolves[1, numVar + 1],
    timeElapsed = Sys.time() - t0
  ))
}
