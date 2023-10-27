### Obsolete constants
NUM_ATLAS_GRAPHS = 1252L
MAX_ATLAS_SIZE   = 7L

### This function parses a written proof in LaTeX format, providing the totals. 
### Note: replace all \ with \\ in the LaTeX source before passing it as input!
### NOTE: This function is no longer used; it is now subsumed by optimiseProofs
parseProof = function(String, N = 9) {
  String %<>%
    str_remove_all("\\{") %>%
    str_remove_all("\\}") %>%
    str_remove_all("\\geq")
  String %<>%
    str_split('\n') %>%
    unlist() %>%
    str_remove_all("X\\[") %>%
    str_remove_all("\\]") %>%
    str_remove_all("[&]") %>%
    str_remove_all("\\\\")
  pairs = String %>%
    str_extract_all("[(][0-9,]+[)]")
  coeffs = String %>%
    str_extract("^[-]?[0-9]+") %>%
    as.integer()
  unitPos = which(is.na(coeffs))
  if (length(unitPos) > 0) {
    coeffs[unitPos] = ifelse(str_starts(String[unitPos], "-"), -1L, 1L)
  }
  rhs = String %>%
    str_extract("[-]?[0-9]+$") %>%
    as.integer()
  Mats = lapply(pairs, function(x) {
    x %<>%
      str_remove_all("[(]") %>%
      str_remove_all("[)]") %>% 
      str_split_fixed(",", n = 2)
    mode(x) = "integer"
    x
  })
  Mat = matrix(0, N, N)
  for (index in 1:length(Mats)) {
    curPairs = Mats[[index]]
    Mat[curPairs] %<>%
      add(coeffs[index])
  }
  totalRhs = sum(rhs)
  sumMat = which(Mat != 0, arr.ind = TRUE) %>%
    as_tibble() %>%
    mutate(val = Mat[cbind(row, col)]) %>%
    arrange(val)
  output = list(sumMat, totalRhs)
  output
}

### This function tests a conjecture that the transitive reduction of the set of
### precedence constraints for the lexicographically smallest element of a coset
### of the automorphism group of a graph is a disjoint union of directed trees. 
### It returns an all-TRUE vector iff the conjecture holds for graphs of order N
### NOTE: This conjecture has now been proven to hold for all permutation groups
testConjecture = function(N) { 
  Graphs = getGraphs(numVerts = N)
  count = length(Graphs)
  fullTest = rep(FALSE, count)
  print(paste("There are", count, "graphs to process"))
  for (ind in 1:count) {
    if (ind %% 10000 == 0) { print(ind) }
    Orbit = getOrbits(Graphs[[ind]])
    fullTest[ind] = !(any(duplicated(Orbit$comps[,2])))
  }
  fullTest
}

### This function constructs the transitive reduction of an input graph
### The graphs' input and output formats are a two-column list of edges
### NOTE: This function is no longer used due to a faster alternative.
getTransitiveReduction = function(edgeList) {
  G0 = graph_from_edgelist(edgeList, directed = TRUE)
  M0 = get.adjacency(G0, sparse = FALSE)
  D0 = distances(G0, mode = "out")
  M1 = (is.finite(D0) & D0 > 0)
  M2 = pmin(M0 %*% M1, 1)
  G1 = graph_from_adjacency_matrix(M2)
  G2 = graph.difference(G0, G1)
  newEdgeList = get.edgelist(G2)
  newEdgeList
}

### This function extracts small graphs from the graph atlas provided in igraph
### Note: it skips graph 0 (empty graph) and graph 1 (the single vertex graph)!
### If extremeOnly = TRUE, only keeps the graphs with sizes within 1 of min/max.
### If symOnly = TRUE, only keeps those with a smaller than average orbit (per order)
getSmallGraphs = function(numGraphs = NUM_ATLAS_GRAPHS, symOnly = FALSE, extremeOnly = FALSE) {
  fn = paste0("AtlasGraphs", ifelse(extremeOnly, "Extreme", ""), ifelse(symOnly, "Sym", ""), ".RData")
  if (!file.exists(fn)) {
    allG = vector("list", numGraphs)
    pos = 1
    for (ind in 1:numGraphs) {
      curG = graph_from_atlas(ind)
      if (is.connected(curG) & gorder(curG) > 1) {
        allG[[pos]] = curG
        pos = pos + 1
      }
    }
    allG = allG[1:(pos - 1)]
    if (extremeOnly) {
      orders = sapply(allG, gorder)
      sizes  = sapply(allG, gsize)
      goodGraphs = which(sizes <= orders | sizes >= choose(orders, 2) - 1)
      allG = allG[goodGraphs]
    }
    if (symOnly) {
      numAutos = as.integer(sapply(allG, automorphisms)["group_size",])
      graphTab = tibble(index = 1:(pos - 1), order = sapply(allG, gorder), orbit = factorial(order)/numAutos)
      finalGraphs = graphTab %>%
        group_by(order) %>%
        mutate(meanOrbit = mean(orbit), betterA = (orbit < meanOrbit)) %>%
        ungroup() %>%
        filter(betterA) %>%
        pull(index)
      allG %<>% magrittr::extract(finalGraphs)
    } 
    save(allG, file = fn)
  }
  else {
    e = new.env()
    load(fn, envir = e)
    allG = get("allG", envir = e)
  }
  allG
}

### Obsolete code from getMediumGraphs
numGraphs = L/4
allG      = vector("list", numGraphs)
for (ind in 1:numGraphs) {
  if (ind %% 10000 == 0) {
    print(ind)
  }
  allG[[ind]] = curText[ind] %>% 
    str_split("  ") %>% 
    unlist %>% 
    str_split_fixed(" ", n = 2) %>% 
    as.numeric() %>% 
    matrix(ncol = 2) %>%
    t() %>%
    magrittr::add(1) %>%
    graph(., directed = FALSE)
}
if (symOnly) {
  numAutos = as.integer(sapply(allG, automorphisms)["group_size",])
  numPerms = factorial(numVerts)
  graphTab = tibble(index = 1:numGraphs, orbit = numPerms/numAutos)
  finalGraphs = graphTab %>%
    mutate(meanOrbit = mean(orbit), betterA = (orbit < meanOrbit)) %>%
    filter(betterA) %>%
    pull(index)
  allG %<>% magrittr::extract(finalGraphs)
}

### Obsolete code from iterateLPProof
allGraphs = getSmallGraphs(symOnly = symOnly, extremeOnly = extremeOnly)
if (maxOrder > MAX_ATLAS_SIZE) {
  for (medSize in (MAX_ATLAS_SIZE + 1):maxOrder) {
    allGraphs %<>% c(getMediumGraphs(medSize, extremeOnly = extremeOnly))
  }
}

graphTab    = tibble(size = sapply(allGraphs, gsize), order = sapply(allGraphs, gorder), n_auto = as.integer(sapply(allGraphs, automorphisms)["group_size",]))

else {
  curSupport   = curSupport + 1
  if (curSupport >= curOrder + 3 && curOrder < maxOrder) {
    curOrder   = curOrder + 1
    curSupport = curOrder
    break
  }
}

### Obsolete code from permuteGraph
### If pairsOnly = TRUE,  returns a list of nC2 adjacency matrices, one per pair
permuteGraph = function(Graph, pairsOnly = FALSE) {
  if (pairsOnly) {
    Graph %<>% get.adjacency(type = "both", sparse = FALSE)
    n = nrow(Graph)
    stopifnot(ncol(Graph) == n)
    allPairs = combn2(1:n)
    allPerms = lapply(1:nrow(allPairs), function(x) { cur = allPairs[x,]; Graph[rev(cur),] = Graph[cur,]; Graph[, rev(cur)] = Graph[, cur]; Graph })
  }
  allPerms
}

tailEdges = allPerms[Graph[, 1], , drop = FALSE]
headEdges = allPerms[Graph[, 2], , drop = FALSE]
swapPos   = which(tailEdges > headEdges)
tempEdges = tailEdges[swapPos]
tailEdges[swapPos] = headEdges[swapPos]
headEdges[swapPos] = tempEdges
allPerms = rbind(tailEdges, headEdges)

### This function optimises the number of red edges in a graph, represented by its edge 
### indices, over the Ramsey polytope on numVertices vertices (so 0 <= curEdges <= nC2)
solveRamseyLP = function(prep, numVertices, curEdges, lower = TRUE, GLPK = FALSE) {
  numVars = choose(numVertices, 2)
  objVec = rep(0, numVars)
  objVec[curEdges] = 1
  if (GLPK) {
    fname = "TempFile.lp"
    numConst = length(prep$dir)
    numNonZeros = tail(prep$mat@p,1)
    model = glpkAPI::initProbGLPK()
    glpkAPI::setProbNameGLPK(model, "Ramsey Polytope Optimisation")
    glpkAPI::setObjDirGLPK(model, ifelse(lower, glpkAPI::GLP_MIN, glpkAPI::GLP_MAX))
    glpkAPI::addColsGLPK(model, ncols = numVars)
    glpkAPI::setColsBndsObjCoefsGLPK(model, j = seq_len(numVars), lb = rep(0, numVars), ub = rep(1, numVars), obj_coef = objVec, type = rep(glpkAPI::GLP_DB, numVars))
    glpkAPI::setColKindGLPK(model, j = seq_len(numVars), kind = rep(glpkAPI::GLP_CV, numVars))
    glpkAPI::addRowsGLPK(model, nrows = numConst)
    glpkAPI::loadMatrixGLPK(model, ne = numNonZeros, ia = prep$mat@i + 1, ja = rep(1:numVars, diff(prep$mat@p)), ra = rep(1, numNonZeros))
    rowTypes = ifelse(prep$dir == "G", glpkAPI::GLP_LO, glpkAPI::GLP_UP)
    glpkAPI::setRowsBndsGLPK(model, i = seq_len(numConst), lb = ifelse(prep$dir == "L", 0, prep$rhs), ub = ifelse(prep$dir == "G", 0, prep$rhs), type = rowTypes)
    glpkAPI::writeLPGLPK(model, fname = fname)
    ### CONTINUE FROM HERE TO CALLING CPLEX AND EXTRACTING A SOLUTION (OBJECTIVE + DUAL VARIABLES)
  }
  solution
}

### This function parses the results of an LP optimisation stored in fname;
### It returns the objective value and the dual vector only at the moment.
### If removeFile = TRUE, the file gets removed after being processed.
parseResults = function(fname, removeFile = FALSE) {
  Lines = readLines(fname) %>%
    str_trim()
  objLine = Lines %>%
    str_detect("objectiveValue") %>%
    which %>%
    min
  objValue = Lines[objLine] %>%
    str_remove("objectiveValue=") %>%
    str_remove_all('\"') %>%
    as.numeric()
  dualLines = Lines %>%
    str_detect("<constraint name") %>%
    which
  dualPos = Lines[dualLines] %>%
    str_extract('\"r_[0-9]+\"') %>%
    str_remove("r_") %>%
    str_remove_all('\"') %>%
    as.numeric()
  stopifnot(all(dualPos == 1:length(dualPos)))
  dualValues = Lines[dualLines] %>%
    str_extract('dual=\"[-]?[0-9]+[\\.]?[0-9]*([eE][-]?[0-9]+)?\"') %>%
    str_remove("dual=") %>%
    str_remove_all('\"') %>%
    as.numeric()
  stopifnot(all(!is.na(dualValues)))
  if (removeFile) {
    file.remove(fname)
  }
  output = list(obj = objValue, duals = dualValues)
  output
}

### Code from permuteGraph: eliminating duplicates obsolete as the group is known
### There are duplicates; radix-sort adjacency lists and eliminate
m = nrow(Graph)
firstOrder  = cbind(as.vector(apply(headEdges, 2, order)), rep(1:numPerms, each = m))
tailEdges   = matrix(tailEdges[firstOrder], ncol = numPerms)
headEdges   = matrix(headEdges[firstOrder], ncol = numPerms)
secondOrder = cbind(as.vector(apply(tailEdges, 2, order)), rep(1:numPerms, each = m))
tailEdges   = matrix(tailEdges[secondOrder], ncol = numPerms)
headEdges   = matrix(headEdges[secondOrder], ncol = numPerms)
allPerms    = rbind(tailEdges, headEdges)
allPerms    = allPerms[, !duplicated(t(allPerms)), drop = FALSE]
stopifnot(ncol(allPerms) == expectedNumPerms)

### Obsolete code from findNextGraph (now called tightenBounds)
model = glpkAPI::initProbGLPK()
glpkAPI::setProbNameGLPK(model, "Ramsey Polytope Optimisation")
glpkAPI::addColsGLPK(model, ncols = numVars)
glpkAPI::setColsBndsGLPK(model, j = seq_len(numVars), lb = rep(0, numVars), ub = rep(1, numVars), type = rep(glpkAPI::GLP_DB, numVars))
glpkAPI::setColKindGLPK(model, j = seq_len(numVars), kind = rep(glpkAPI::GLP_CV, numVars))
glpkAPI::addRowsGLPK(model, nrows = numConst)
glpkAPI::loadMatrixGLPK(model, ne = numNonZeros, ia = prep$mat@i + 1, ja = rep(1:numVars, diff(prep$mat@p)), ra = rep(1, numNonZeros))
rowTypes = ifelse(prep$dir == "G", glpkAPI::GLP_LO, glpkAPI::GLP_UP)
glpkAPI::setRowsBndsGLPK(model, i = seq_len(numConst), lb = ifelse(prep$dir == "L", 0, prep$rhs), ub = ifelse(prep$dir == "G", 0, prep$rhs), type = rowTypes)
glpkAPI::writeLPGLPK(model, fname = fname)
### Parsing section
curResultL = parseResults(curOutFiles[1])
curObjL    = curResultL$obj
boundL     = as.integer(ifelse(near(curObjL, round(curObjL)), round(curObjL), ceiling(curObjL)))
curDualsL  = curResultL$duals ## curResultL$extra$lambda
proofL     = NULL
if (boundL >= 1               && !(all(near(curDualsL[1:numConst], 0)))) {
  proofL   = constructProof(prep, curDualsL[1:numConst], lower = TRUE, short = short)
}
if (symRS) {
  boundU = curSize - boundL
  proofU = proofL
} else {
  curResultU = parseResults(curOutFiles[2])
  curObjU    = curResultU$obj
  boundU     = as.integer(ifelse(near(curObjU, round(curObjU)), round(curObjU), floor(curObjU)))
  curDualsU  = curResultU$duals ## curResultU$extra$lambda
  proofU     = NULL
  if (boundU <= curRow$size - 1 && !(all(near(curDualsU[1:numConst], 0)))) {
    proofU   = constructProof(prep, curDualsU[1:numConst], lower = FALSE, short = short)
  }
}

### This function finds bounds on the number of red edges in a graph; NOTE: has been renamed to tightenBounds in Ramsey.R!
### The Ramsey polytope is constructed over the specified order, numVertices; known bounds are provided via boundTab
### The graph will be the first (respectively, all) of the input graphs provided (allGraphs) to yield a non-trivial bound
### The allGraphAdj variable contains all the distinct permuted versions of each graph's adjacency list in the same order
### that also contain relevant information in graphInfo (in particular, whose indices are contained in its first column).
findNextGraph = function(numVertices, allGraphs, allGraphAdj, boundTab, graphInfo) {
  L = nrow(graphInfo)
  prep        = prepareRamseyLP(numVertices = numVertices, allGraphAdj = allGraphAdj, graphBounds = boundTab, graphInfo = graphInfo, sparse = TRUE)
  numConst    = nrow(prep$mat)
  numVars     = choose(numVertices, 2)
  print(paste("There are", numConst, "constraints over", numVars, "variables"))
  output      = vector("list", L)
  auxMatrix   = createPositionMatrix(numVertices)
  print(paste("There are", L, "graphs to process"))
  control    = list(trace = 0, preind = 0, method = 2)
  objVec     = rep(0, numVars)
  for (ind in 1:L) {
    print(ind)
    curResult  = NULL
    curRow     = graphInfo %>% slice(ind)
    curGraph   = allGraphs[[curRow$index]]
    curEdges   = auxMatrix[as_edgelist(curGraph)]
    curObjVec  = objVec
    curObjVec[curEdges] = 1
    curResultL = Rcplex(cvec = curObjVec, Amat = prep$mat, bvec = prep$rhs, lb = 0, ub = 1, objsense = "min", sense = prep$dir, control = control)
    curObjL    = curResultL$obj
    boundL     = as.integer(ifelse(near(curObjL, round(curObjL)), round(curObjL), ceiling(curObjL)))
    curDualsL  = curResultL$extra$lambda
    proofL     = NULL
    if (boundL >= 1               && !(all(near(curDualsL[1:numConst], 0)))) {
      proofL   = constructProof(prep, curDualsL[1:numConst], lower = TRUE, short = TRUE)
    }
    curResultU = Rcplex(cvec = curObjVec, Amat = prep$mat, bvec = prep$rhs, lb = 0, ub = 1, objsense = "max", sense = prep$dir, control = control)
    curObjU    = curResultU$obj
    boundU     = as.integer(ifelse(near(curObjU, round(curObjU)), round(curObjU), floor(curObjU)))
    curDualsU  = curResultU$extra$lambda
    proofU     = NULL
    if (boundU <= curRow$size - 1 && !(all(near(curDualsU[1:numConst], 0)))) {
      proofU   = constructProof(prep, curDualsU[1:numConst], lower = FALSE, short = TRUE)
    }
    if (!is.null(proofL) || !is.null(proofU)) {
      curResult = list(graph = curRow$index, proofL = proofL, proofU = proofU, boundL = boundL, boundU = boundU)
    }
    output[[ind]] = curResult
  }
  output
}

### This function optimises the number of red edges in a graph, represented by its edge 
### indices, over the Ramsey polytope on numVertices vertices (so 0 <= curEdges <= nC2)
solveRamseyLP = function(prep, numVertices, curEdges, lower = TRUE, GLPK = FALSE) {
  numVars = choose(numVertices, 2)
  objVec = rep(0, numVars)
  objVec[curEdges] = 1
  if (GLPK) {
    fname = "TempFile.lp"
    numConst = length(prep$dir)
    numNonZeros = tail(prep$mat@p,1)
    model = glpkAPI::initProbGLPK()
    glpkAPI::setProbNameGLPK(model, "Ramsey Polytope Optimisation")
    glpkAPI::setObjDirGLPK(model, ifelse(lower, glpkAPI::GLP_MIN, glpkAPI::GLP_MAX))
    glpkAPI::addColsGLPK(model, ncols = numVars)
    glpkAPI::setColsBndsObjCoefsGLPK(model, j = seq_len(numVars), lb = rep(0, numVars), ub = rep(1, numVars), obj_coef = objVec, type = rep(glpkAPI::GLP_DB, numVars))
    glpkAPI::setColKindGLPK(model, j = seq_len(numVars), kind = rep(glpkAPI::GLP_CV, numVars))
    glpkAPI::addRowsGLPK(model, nrows = numConst)
    glpkAPI::loadMatrixGLPK(model, ne = numNonZeros, ia = prep$mat@i + 1, ja = rep(1:numVars, diff(prep$mat@p)), ra = rep(1, numNonZeros))
    rowTypes = ifelse(prep$dir == "G", glpkAPI::GLP_LO, glpkAPI::GLP_UP)
    glpkAPI::setRowsBndsGLPK(model, i = seq_len(numConst), lb = ifelse(prep$dir == "L", 0, prep$rhs), ub = ifelse(prep$dir == "G", 0, prep$rhs), type = rowTypes)
    glpkAPI::writeLPGLPK(model, fname = fname)
    ### CONTINUE FROM HERE TO CALLING CPLEX AND EXTRACTING A SOLUTION (OBJECTIVE + DUAL VARIABLES)
  } else {
    control = list(trace = 0, preind = 0, method = 2)
    solution = Rcplex(cvec = objVec, Amat = prep$mat, bvec = prep$rhs, lb = 0, ub = 1, objsense = ifelse(lower, "min", "max"), sense = prep$dir, control = control)
  }
  solution
}

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

### Obsolete functionality used to construct LPs via GLPK
outputFiles  = paste0("Solution", nRound, "G", outer(c("Lower", "Upper"), 1:L, function(x,y) {paste0(y,x)}), ".sol")
startLine    = 'source ~/.bash_profile; ./cplex -c'
settingLines = c(paste('re', fname), paste('set', 'prep', 'pres', 'n', collapse = ' '), paste('set', 'lpm', 2, collapse = ' '))
finalLine    = 'qu'
initDir = getwd()
setwd(CPLEX_DIR)
if (!file.exists(curOutFiles[1])) {
  changeLines  = rep("", curSize + ifelse(symRS, 2, 5))
  for (pos in 1:curSize) {
    changeLines[pos] = paste('ch', 'ob', paste0('X', curEdges[pos]), 1, collapse = ' ')
  }
  changeLines[curSize + (1:2)] = c('op', paste('wr', curOutFiles[1]))
  ## if (!symRS) {
  changeLines[curSize + (3:5)] = c(paste('ch', 'se', '0', 'max', collapse = ' '), 'op', paste('wr', curOutFiles[2]))
  ## }
  fullCommand  = system(paste(startLine, paste(map_chr(c(settingLines, changeLines, finalLine), ~{paste('"', ., '"', collapse = "")}), collapse = " ")))
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

### This obsolete function parses the optimal value of an LP optimisation stored in fname
parseObjective = function(fname) {
  myResult = system(paste('head', '-5', fname, collapse = ' '), intern = TRUE)
  objValue = myResult[5] %>%
    str_remove("objectiveValue=") %>%
    str_remove_all('\"') %>%
    as.numeric()
  objValue
}

### This obsolete function parses the dual values of an LP optimisation stored in fname
parseDuals = function(fname) {
  Lines = readLines(fname) %>%
    str_trim()
  dualLines = Lines %>%
    str_detect("^<constraint name") %>%
    which
  dualValues = Lines[dualLines] %>%
    str_extract('dual=\"[-]?[0-9]+[\\.]?[0-9]*([eE][-]?[0-9]+)?\"') %>%
    str_remove("dual=") %>%
    str_remove_all('\"') %>%
    as.numeric()
  stopifnot(all(!is.na(dualValues)))
  dualValues
}
