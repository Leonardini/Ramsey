### This code used to be part of the prepareRamseyLP function
if (curLength == 2) {
  allPermEdges = matrix(allPermEdges, ncol = curNOpts)
  allPermEdges = as.vector(rbind(allPermEdges, allPermEdges))
}
print("Sorting the position matrix by column")
myOrder = order(Mat[, 2])
print("Creating the triplet matrix")
Mat = simple_triplet_matrix(i = Mat[myOrder, 1], j = Mat[myOrder, 2], v = rep(1, length(myOrder)), nrow = numConstraints, ncol = numVars)
### Alternative using the sparseMatrix functionality instead
myOrder = sort.int(Mat[, 2], method = "radix", index.return = TRUE)
Mat = sparseMatrix(i = Mat[myOrder$ix, 1], j = myOrder$x, x = 1, dims = c(numConstraints, numVars), index1 = TRUE)

### This code used to be part of the identifySubsumedBounds function
subIso = matrix(FALSE, L, L)
uniqueGraphs = fullTab$graph
testGraphs   = allGraphs[uniqueGraphs]
pos = 1
while (pos < L) {
  curGraph = testGraphs[[1]]
  testGraphs %<>% magrittr::extract(-1)
  subIso[pos, (pos + 1):L] = sapply(testGraphs, function(x) { subgraph_isomorphic(curGraph, x) })
  pos = pos + 1
}
curSuperIso = uniqueGraphs[subIso[index, ]]
curSubIso   = uniqueGraphs[subIso[, index]]
subsumedBySuper = fullTab %>%
  filter(graph %in% curSuperIso & (bound - curRow$bound) * ifelse(lower, 1, -1) >= ifelse(rep(lower, K), size - curRow$size, rep(0, K)))
subsumedBySub   = fullTab %>%
  filter(graph %in% curSubIso   & (bound - curRow$bound) * ifelse(lower, 1, -1) >= ifelse(rep(lower, K), rep(0, K), curRow$size - size))
fullTab$subsumed[index] = (nrow(subsumedBySuper) > 0 || nrow(subsumedBySub) > 0)

### This function computed a vector of edge positions corresponding to a graph
### It has been superseded by a more direct computation using graph matrices
### The edges are ordered as follows: 12, 13, ..., 1N, 23, ..., 2N, ..., (N-1)N,
### where N is the integer specified by numVertices
### If Map is provided, it is an increasing vector of integers, one per vertex
### If so, the vertices are first mapped according to it before naming the edges
getEdgePositions = function(Graph, numVertices, Map = NULL) {
  offsets = c(0, cumsum(numVertices - (1:(numVertices - 1))))
  if ("igraph" %in% class(Graph)) {
    Graph %<>%
      as_adjacency_matrix(type = "upper", sparse = FALSE)
  }
  allEdges = which(Graph == 1, arr.ind = TRUE) %>%
    as_tibble() %>%
    filter(row < col) %>%
    arrange(row, col)
  if (!is.null(Map)) {
    stopifnot(all(Map <= numVertices))
    allEdges = allEdges %>% 
      mutate_all(~{ Map[.] })
  }
  allEdges %<>%
    mutate(delta = col - row, pos = offsets[row] + delta) %>%
    pull(pos)
  allEdges
}

### This function prepared to optimise over the Ramsey polytope with specified bounds
### Unfortunately, it had issues with not working correctly when graphs were repeated
prepareRamseyLP = function(numVertices, allGraphAdj, graphBounds, graphInfo, sparse = TRUE) {
  numVars = choose(numVertices, 2)
  allVars = 1:numVars
  graphBounds %<>%
    filter(!subsumed) %>%
    inner_join(graphInfo, by = join_by(graph == index)) %>%
    mutate(numOpts = choose(numVertices, order), numCombinations = n_auto * numOpts, fullCombinations = numCombinations * size)
  lowerBounds = graphBounds %>%
    filter(direction == LOWER)
  upperBounds = graphBounds %>%
    filter(direction == UPPER)
  numCombinationsL = lowerBounds$numCombinations
  numCombinationsU = upperBounds$numCombinations
  numConstraintsL  = sum(numCombinationsL)
  numConstraintsU  = sum(numCombinationsU)
  Dir = c(rep("G", numConstraintsL), rep("L", numConstraintsU))
  Rhs = c(rep(lowerBounds$bound, numCombinationsL), rep(upperBounds$bound, numCombinationsU))
  if (sparse) {
    matLower = matrix(0, sum(lowerBounds$fullCombinations), 2)
    matUpper = matrix(0, sum(upperBounds$fullCombinations), 2)
  } else {
    matLower = matrix(0, numConstraintsL, numVars, dimnames = list(c(), allVars))
    matUpper = matrix(0, numConstraintsU, numVars, dimnames = list(c(), allVars))
  }
  posL = 1
  posU = 1
  cntL = 0
  cntU = 0
  MapL = cbind(1:numConstraintsL, 0)
  MapU = cbind(1:numConstraintsU, 0)
  for (ind in unique(graphBounds$graph)) {
    curGraphs = allGraphAdj[[ind]]
    curInfo   = graphBounds %>%
      filter(graph == ind)
    curOrder  = curInfo$order[1]
    curSize   = curInfo$size[1]
    curOpts   = combn(1:numVertices, curOrder)
    if (curOrder == numVertices) {
      curOpts = matrix(curOpts, ncol = 1)
    }
    for (index in 1:length(curGraphs)) {
      myGraph = curGraphs[[index]]
      for (j in 1:ncol(curOpts)) {
        curEdges = getEdgePositions(myGraph, numVertices = numVertices, Map = curOpts[,j])
        if (LOWER %in% curInfo$direction) {
          curN = curInfo %>%
            filter(direction == LOWER) %>%
            pull(number)
          if (sparse) {
            matLower[cntL + (1:curSize), ] = cbind(rep(posL, curSize), curEdges)
            cntL = cntL + curSize
          } else {
            MatL[posL, curEdges] = 1
          }
          MapL[posL, 2] = curN
          posL = posL + 1
        }
        if (UPPER %in% curInfo$direction) {
          curN = curInfo %>%
            filter(direction == UPPER) %>%
            pull(number)
          if (sparse) {
            matUpper[cntU + (1:curSize), ] = cbind(rep(posU, curSize), curEdges)
            cntU = cntU + curSize
          } else {
            matUpper[posU, curEdges] = 1
          }
          MapU[posU, 2] = curN
          posU = posU + 1
        }
      }
    }
  }
  if (sparse) {
    matUpper[, 1] %<>% 
      add(numConstraintsL)
    fullMat = rbind(matLower, matUpper)
    myOrder = order(fullMat[, 2])
    fullMat = fullMat[myOrder, ]
    Mat = simple_triplet_matrix(i = fullMat[, 1], j = fullMat[, 2], v = rep(1, length(myOrder)), nrow = numConstraintsL + numConstraintsU, ncol = numVars)
  } else {
    Mat = rbind(matLower, matUpper)
  }
  MapU[, 1] %<>%
    add(numConstraintsL)
  Map = rbind(MapL, MapU)
  output = list(mat = Mat, dir = Dir, rhs = Rhs, map = Map)
  output
}

### Obsolete code from identifySubsumedBounds (because both subgraph and supergraph isomorphisms are relevant)
for (pos in 1:L) {
  curTest = fullTab %>%
    slice(pos)
  curGraph = allGraphs[[curTest$graph]]
  if (lower) {
    curCandidates = fullTab %>%
      filter(size <= curTest$size & order <= curTest$order & bound <= curTest$bound & graph != curTest$graph)
  } else {
    curCandidates = fullTab %>%
      filter(size >= curTest$size & order >= curTest$order & bound >= curTest$bound & graph != curTest$graph)
  }
  if (nrow(curCandidates) > 0) {
    curCandGraphs = allGraphs[curCandidates$graph]
    if (lower) {
      curCandidates %<>%
        mutate(super_iso = sapply(curCandGraphs, function(x) {subgraph_isomorphic(x, curGraph)})) %>%
        mutate(subsuming = super_iso   & (bound - curTest$bound) >= (size - curTest$size))
    } else {
      curCandidates %<>%
        mutate(sub_iso = sapply(curCandGraphs, function(y) {subgraph_isomorphic(curGraph, y)})) %>%
        mutate(subsuming = sub_iso     & (curTest$bound - bound) >= (curTest$size - size))
    }
    if (any(curCandidates$subsuming)) {
      fullTab$subsumed[pos] = TRUE
    }
  }
}

### Omitted functionality from iterateLPProof; the proofs are now not checked due to numerical issues
### Former function comment: If checkProofs is TRUE, independently checks each new constraint's accuracy
if (checkProofs) {
  curEdges   = getEdgePositions(curResult$graph, numVertices = curSupport, Map = NULL)
  proofSumL  = colSums(curProofL * curProofL[,1])[-1]
  posSumL    = sum((proofSumL * ifelse(proofSumL > 0, 1, 0))[-curEdges]) - proofSumL["rhs"]
  proofLower = computeBestBound(proofSumL[curEdges], proofSumL["rhs"] - posSumL, lower = TRUE)
  if (!near(proofLower, lowerBound)) {
    print(curProofL)
    print(curEdges)
    print(proofSumL)
    print(lowerBound)
    stop("Erroneous lower bound!")
  }
}

if (checkProofs) {
  curEdges   = getEdgePositions(curResult$graph, numVertices = curSupport, Map = NULL)
  proofSumU  = colSums(curProofU * curProofU[,1])[-1]
  negSumU    = sum((proofSumU * ifelse(proofSumU < 0, -1, 0))[-curEdges])
  proofUpper = computeBestBound(proofSumU[curEdges], proofSumU["rhs"] + negSumU, lower = FALSE)
  if (!near(proofUpper, upperBound)) {
    print(curProofU)
    print(curEdges)
    print(proofSumU)
    print(upperBound)
    stop("Erroneous upper bound!")
  }
}


### This function computes the best (lower if lower = TRUE, upper otherwise) bound 
### for the sum of N [0,1] variables from given coefficient vector and constant value
### For instance, given 3x + 4y + 5z <= 6 with x, y, z in [0, 1], the best bound can
### be obtained by adding x <= 1 to 3x + 4y + 4z <= 6 to get 7/4 = 1.75, rounded to 2
computeBestBound = function(coeffVector, constValue, lower = TRUE) {
  if (!(all(near(coeffVector, round(coeffVector))))) {
    print("Warning: the vector is not close to an integer vector!")
    print(coeffVector)
    print(constValue)
  } else {
    coeffVector %<>%
      as.integer()
  }
  RLE = coeffVector %>%
    sort(decreasing = lower) %>%
    rle()
  mults    = RLE$lengths
  values   = RLE$values
  options  = (constValue - cumsum(mults * values))/values + cumsum(mults)
  if (lower) {
    output = as.integer(ceiling(max(options)))
  } else {
    output = as.integer(floor  (min(options)))
  }
  output
}

### This function computes a vector of edge names corresponding to a graph
### If Map is provided, it is an increasing vector of integers, one per vertex
### If so, the vertices are first mapped according to it before naming the edges
getEdgeNames = function(Graph, Map = NULL) {
  if (class(Graph)[1] == "igraph") {
    Graph %<>%
      get.adjacency(type = "upper", sparse = FALSE)
  }
  allEdges = which(Graph == 1, arr.ind = TRUE) %>%
    as_tibble() %>%
    filter(row < col) %>%
    arrange(row, col)
  if (!is.null(Map)) {
    allEdges = allEdges %>% 
      mutate_at(c("row", "col"), ~{ Map[.] })
  }
  allNames = paste0(allEdges$row, "C", allEdges$col)
  allNames
}

### Obsolete function for processing the Ramsey LP/ILP via the GLPK interface 
processInputFiles = function(modelFile = "RamseyModel", dataFile = "RamseyData", inputFile = NULL, LPonly = TRUE) {
  problem <- glpkAPI::initProbGLPK()
  pointer <- glpkAPI::mplAllocWkspGLPK(ptrtype = "tr_wksp")
  glpkAPI::mplReadModelGLPK(wk = pointer, fname = modelFile, skip = 1)
  if (!is.null(dataFile)) {
    print("Processing the data in the specified input file")
    glpkAPI::mplReadDataGLPK(wk = pointer, fname = dataFile)
  }
  glpkAPI::mplGenerateGLPK(wk = pointer)
  glpkAPI::mplBuildProbGLPK(wk = pointer, lp = problem)
  glpkAPI::setMIPParmGLPK(glpkAPI::PRESOLVE, glpkAPI::GLP_ON)
  glpkAPI::setMIPParmGLPK(glpkAPI::MSG_LEV, glpkAPI::GLP_MSG_ALL)
  numCol <- glpkAPI::getNumColsGLPK(problem)
  varNames <- rep("", numCol)
  for (j in 1:numCol) {varNames[j] <- glpkAPI::getColNameGLPK(problem, j)}
  if (!is.null(inputFile)) {
    print("Saving the input to an LP file")
    glpkAPI::writeLPGLPK(problem, fname = inputFile)
  }
  if (LPonly) {
    outputFN <- gsub(".lp", ".sol", inputFile)
    if (!outputFN %in% list.files()) {
      runTime <- system.time(system2("cbc", args = c("-import", inputFile, "-solve", "-solution", outputFN)))
      fullSolution <- CBC_parse(outputFN)
      solution <- fullSolution[[1]]
      objValue <- fullSolution[[2]]
    }
  } else {
    runTime  <- system.time({capture.output(glpkAPI::solveMIPGLPK(problem))})
    solution <- glpkAPI::mipColsValGLPK(problem)
    names(solution) <- varNames
    solution %<>% 
      postprocessSolutionx
    objValue <- glpkAPI::mipObjValGLPK(problem)
    glpkAPI::mplFreeWkspGLPK(pointer)
    glpkAPI::delProbGLPK(problem)
  }
  output <- list(time = runTime, solution = solution, objective = objValue)
  output
}

### obsolete parts of the code dealing with the upper bounds
if (!is.null(curProofU)) {
  proofSumU  = colSums(curProofU * curProofU[,1])[-1]
  proofSumU[near(proofSumU, 0)] = Inf
  proofMinU  = min(proofSumU[curEdges])
  proofUpper = floor(proofSumU["rhs"]/proofMinU)
  stopifnot(near(proofUpper, upperBound))
}

curResultL = curResults$lower
curResultU = curResults$upper
curObjL = curResultL$objval
curObjU = curResultU$objval
proofL = NULL
proofU = NULL
boundL = ifelse(near(curObjL, round(curObjL)), as.integer(round(curObjL)), as.integer(ceiling(curObjL)))
boundU = ifelse(near(curObjU, round(curObjU)), as.integer(round(curObjU)), as.integer(floor  (curObjU)))
if (boundL >= 1        && sum(!near(curResultL$duals[1:numConst], 0)) > 1) { ### ignore results due to subgraphs
  print(paste("Found a lower bound of", boundL, "at index", ind))
  proofL = constructProof(prep$mat, prep$dir, prep$rhs, curResultL$duals, lower = TRUE)
} else {
  boundL = -Inf
}
if (boundU <= curM - 1 && sum(!near(curResultU$duals[1:numConst], 0)) > 1) { ### ignore results due to subgraphs
  print(paste("Found an upper bound of", boundU, "at index", ind))
  proofU = constructProof(prep$mat, prep$dir, prep$rhs, curResultU$duals, lower = FALSE)
} else {
  boundU = Inf
}
if (!is.null(proofL) || !is.null(proofU)) {
  curResult = list(ind = ind, graph = curGraph, proofU = proofU, proofL = proofL, boundU = boundU, boundL = boundL)
  if (stopAfterOne) {
    output = list(curResult)
    break
  } else {
    output = c(output, list(curResult))
  }
}

keepRows = which(is.finite(Rhs))
Mat = Mat[keepRows, , drop = TRUE]
Dir = Dir[keepRows]
Rhs = Rhs[keepRows]

SolU = lp("max", objective.in = objVec, const.mat = Mat, const.dir = Dir, const.rhs = Rhs, compute.sens = TRUE)
output = list(lower = SolL, upper = SolU)
output