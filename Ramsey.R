library(combinat)
library(igraph)
library(Rcplex)
library(magrittr)
library(Matrix)
library(slam)
library(tidyverse)

NUM_ATLAS_GRAPHS = 1252L
MAX_ATLAS_SIZE   = 7L

setwd("/Users/lchindelevitch/Downloads/NonPriority/Conjectures/Ramsey numbers")

### TODO: Consider implementing a version that looks for cyclic orientations too
### TODO: The idea is to only constrain delta vectors which admit a bipartition!

### This function extracts small graphs from the graph atlas provided in igraph
### Note: it skips graph 0 (empty graph) and graph 1 (the single vertex graph)!
getSmallGraphs = function(numGraphs = NUM_ATLAS_GRAPHS, symOnly = FALSE) {
  fn = paste0("AtlasGraphs", ifelse(symOnly, "Sym", ""), ".RData")
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
    if (symOnly) {
      numAutos = as.integer(sapply(allG, automorphisms)["group_size",])
      # allPairs = allG %>%
      #   lapply(permuteGraph, pairsOnly = TRUE)
      # numPairs = sapply(1:length(allPairs), function(x) { curG = as_adjacency_matrix(allG[[x]], sparse = FALSE); sum(sapply(allPairs[[x]], function(y) {all(curG == y) })) })
      graphTab = tibble(index = 1:(pos - 1), order = sapply(allG, gorder), orbit = factorial(order)/numAutos) #, twins = numPairs)
      finalGraphs = graphTab %>%
        group_by(order) %>%
        mutate(meanOrbit = mean(orbit), betterA = (orbit < meanOrbit)) %>% # , meanTwins = weighted.mean(twins, orbit), betterB = (twins > meanTwins)) %>%
        ungroup() %>%
        filter(betterA) %>% # & betterB) %>%
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

### This function extracts small graphs from the graph file provided by B. McKay
### Make sure that the conversion is done with the -l0 option in showg_mac64 (!)
getMediumGraphs = function(numVerts = 8, symOnly = FALSE) {
  fn = paste0("Graphs", ifelse(symOnly, "Sym", ""), numVerts, ".RData")
  if (!file.exists(fn)) {
    curText = readLines(paste0("graph", numVerts, "c.txt"))
    L = length(curText)
    goodLines = 4 * (1:(L/4))
    curText = curText[goodLines]
    numGraphs = L/4
    allG      = vector("list", numGraphs)
    for (ind in 1:numGraphs) {
      if (ind %% 1000 == 0) {
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
      # allPairs = allG %>%
      #   lapply(permuteGraph, pairsOnly = TRUE)
      # numPairs = sapply(1:length(allPairs), function(x) { curG = as_adjacency_matrix(allG[[x]], sparse = FALSE); sum(sapply(allPairs[[x]], function(y) {all(curG == y) })) })
      numPerms = factorial(numVerts)
      graphTab = tibble(index = 1:numGraphs, orbit = numPerms/numAutos) # , twins = numPairs)
      finalGraphs = graphTab %>%
        mutate(meanOrbit = mean(orbit), betterA = (orbit < meanOrbit)) %>% # , meanTwins = weighted.mean(twins, orbit), betterB = (twins > meanTwins)) %>%
        filter(betterA) %>% # & betterB) %>%
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

### This function tries to iteratively construct a proof that Ramsey(r,s) <= n
### At the moment, the checks for each order and stops when numVertices >= order + 2
iterateLPProof = function(n = 6L, r = 3L, s = r, symOnly = TRUE, maxOrder = 2 * n / 3, minAuto = factorial(min(r,s))) {
  allGraphs = getSmallGraphs(symOnly = symOnly)
  if (maxOrder > MAX_ATLAS_SIZE) {
    for (medSize in (MAX_ATLAS_SIZE + 1):maxOrder) {
      allGraphs %<>% c(getMediumGraphs(medSize, symOnly = symOnly))
    }
  }
  graphTab    = tibble(size = sapply(allGraphs, gsize), order = sapply(allGraphs, gorder), n_auto = as.integer(sapply(allGraphs, automorphisms)["group_size",]))
  if (!is.null(minAuto)) {
    goodGraphs = which(graphTab$order <= maxOrder & graphTab$n_auto >= minAuto)
  } else {
    goodGraphs  = which(graphTab$order <= maxOrder)
  }
  allGraphs %<>% magrittr::extract(goodGraphs)
  graphTab  %<>% slice(goodGraphs) %>%
    mutate(index = 1:length(allGraphs), .before = 1)
  numGraphs = length(goodGraphs)
  allGraphAdj = vector("list", numGraphs)
  print(paste("There are", numGraphs, "graphs to preprocess!"))
  for (index in 1:numGraphs) {
    allGraphAdj[[index]] = permuteGraph(allGraphs[[index]], pairsOnly = FALSE)
    if (index %% 100 == 0) { print(index) }
  }
  graphTab %<>% mutate(orbit = sapply(allGraphAdj, ncol))
  rInd = graphTab %>%
    filter(order == r & size == choose(r, 2)) %>%
    pull(index)
  sInd = graphTab %>%
    filter(order == s & size == choose(s, 2)) %>%
    pull(index)
  boundTab = tibble(number = 1:2, support = c(r, s), graph = c(rInd, sInd), direction = c("G", "L"), bound = c(1, choose(s, 2) - 1), round = 0L, subsumed = FALSE)
  boundIndex = 3
  proofs     = vector("list", n^2)
  curSupport = min(r, s) + 1
  curOrder   = min(r, s) + 1
  curRound   = 0
  while (curOrder <= maxOrder) {
    curGraphTab = graphTab %>% 
      filter(order >= min(r, s) & order <= curOrder)
    print(paste("Currently looking at graphs with", curOrder, "vertices"))
    while (curSupport <= n) {
      print(paste("Exploring graphs on up to", curOrder, "vertices with support size", curSupport))
      lowerBounds = identifySubsumedBounds(allGraphs, boundTab, curGraphTab, lower = TRUE)
      upperBounds = identifySubsumedBounds(allGraphs, boundTab, curGraphTab, lower = FALSE)
      boundTab    = bind_rows(lowerBounds, upperBounds)
      curRound    = curRound + 1
      fullResult  = findNextGraph(numVertices = curSupport, allGraphs = allGraphs, allGraphAdj = allGraphAdj, boundTab = boundTab, graphInfo = curGraphTab)
      improved   = FALSE
      contradict = FALSE
      if (!all(sapply(fullResult, is.null))) {
        for (ind in 1:length(fullResult)) {
          curResult = fullResult[[ind]]
          if (!is.null(curResult)) {
            curProofL  = curResult$proofL
            curProofU  = curResult$proofU
            lowerBound = curResult$boundL
            upperBound = curResult$boundU
            curGraph   = curResult$graph
            prevLowerBound = max(c(0,                                                        boundTab %>% filter(graph == curGraph & direction == "G") %>% pull(bound)))
            if (lowerBound > prevLowerBound) {
              proofs[[boundIndex]] = curProofL
              boundTab %<>% 
                bind_rows(tibble(number = boundIndex, support = curSupport, graph = curGraph, direction = "G", bound = lowerBound, round = curRound, subsumed = (nrow(curProofL) == 1)))
              boundIndex = boundIndex + 1
            }
            prevUpperBound = min(c(curGraphTab %>% filter(index == curGraph) %>% pull(size), boundTab %>% filter(graph == curGraph & direction == "L") %>% pull(bound)))
            if (upperBound < prevUpperBound) {
              proofs[[boundIndex]] = curProofU
              boundTab %<>% 
                bind_rows(tibble(number = boundIndex, support = curSupport, graph = curGraph, direction = "L", bound = upperBound, round = curRound, subsumed = (nrow(curProofU) == 1)))
              boundIndex = boundIndex + 1
            }
            if (lowerBound > prevLowerBound | upperBound < prevUpperBound) {
              improved = TRUE
              if (lowerBound > upperBound) {
                contradict = TRUE
              }
            }
          }
        }
      }
      if (contradict) {
        print("Found a contradiction!")
        break
      }
      if (!improved) {
        curSupport   = curSupport + 1
        if (curSupport >= curOrder + 3 && curOrder < maxOrder) {
          curOrder   = curOrder + 1
          curSupport = curOrder
          break
        }
      }
    }
    if (contradict) {
      break
    }
  }
  boundTab %<>%
    mutate(across(!direction, as.integer)) %>%
    group_by(round) %>%
    mutate(iter = cur_group_id()) %>%
    ungroup() %>%
    select(-round, -subsumed) %>%
    arrange(number)
  proofs = proofs[1:nrow(boundTab)]
  contradictions = NULL
  if (contradict) {
    contradictions = reconstructContradictions(proofs, boundTab)
  }
  firstGraphs = lapply(allGraphAdj, function(x) { matrix(x[,1], ncol = 2) })
  output = list(proofs = proofs, boundTab = boundTab, graphs = firstGraphs, contradictions = contradictions)
  output
}

### This function finds bounds on the number of red edges in a graph
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
  for (ind in 1:L) {
    print(ind)
    curResult  = NULL
    curRow     = graphInfo %>% slice(ind)
    curGraph   = allGraphs[[curRow$index]]
    curEdges   = auxMatrix[as_edgelist(curGraph)]
    curResultL = solveRamseyLP(prep, numVertices, curEdges, lower = TRUE)
    curObjL    = curResultL$obj
    boundL     = as.integer(ifelse(near(curObjL, round(curObjL)), round(curObjL), ceiling(curObjL)))
    curDualsL  = curResultL$extra$lambda
    proofL     = NULL
    if (boundL >= 1               && !(all(near(curDualsL[1:numConst], 0)))) {
      proofL   = constructProof(prep, curDualsL[1:numConst], lower = TRUE, short = TRUE)
    }
    curResultU = solveRamseyLP(prep, numVertices, curEdges, lower = FALSE)
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

### This function prepares to optimise over the Ramsey polytope with specified bounds
prepareRamseyLP = function(numVertices, allGraphAdj, graphBounds, graphInfo, sparse = TRUE) {
  numVars = choose(numVertices, 2)
  auxMatrix = createPositionMatrix(numVertices = numVertices)
  graphBounds %<>%
    filter(!subsumed) %>%
    inner_join(graphInfo, by = join_by(graph == index)) %>%
    arrange(graph, direction) %>%
    mutate(numOpts = choose(numVertices, order), numCombinations = orbit * numOpts, fullCombinations = numCombinations * size)
  numConstraints = sum(graphBounds$numCombinations)
  Dir = rep("", numConstraints)
  Rhs = rep(0,  numConstraints)
  Map = matrix(0,  numConstraints, 3)
  Mat = matrix(0, ifelse(sparse, sum(graphBounds$fullCombinations), numConstraints), ifelse(sparse, 2, numVars))
  pos = 0
  cnt = 0
  print(paste("There are", length(unique(graphBounds$graph)), "graphs to prepare for the LP"))
  for (ind in unique(graphBounds$graph)) {
    curGraphs = allGraphAdj[[ind]]
    curInfo   = graphBounds %>%
      filter(graph == ind)
    curLength = nrow(curInfo)
    stopifnot(curLength <= 2)
    curSize   = curInfo$size[1]
    curCombos = curInfo$numCombinations[1]
    curNOpts  = curInfo$numOpts[1]
    curNOrbit = curInfo$orbit[1]
    curNConst = curLength * curNOpts
    curNEntry = curNConst * curSize
    curNFull  = curLength * curCombos
    Rhs[pos + (1:curNFull)]   = rep(curInfo$bound,        curCombos)
    Dir[pos + (1:curNFull)]   = rep(curInfo$direction,    curCombos)
    miniIndices = rep(1:curNOpts, each = curLength)
    Map[pos + (1:curNFull), ] = cbind(rep(curInfo$number, curCombos), rep(1:curNOrbit, each = curNConst), rep(miniIndices, curNOrbit))
    curOpts   = combn(1:numVertices, curInfo$order[1])
    if (curNOpts == 1) {
      curOpts = matrix(curOpts, ncol = 1)
    }
    if (curLength == 2) {
      curOpts = curOpts[, miniIndices]
    }
    for (index in 1:ncol(curGraphs)) { ### This updated version processes all the combinations at once!
      curBaseEdges = curGraphs[, index]
      curPermEdges = curOpts[curBaseEdges, , drop = FALSE]
      allPermEdges = auxMatrix[cbind(as.vector(curPermEdges[1:curSize, ]), as.vector(curPermEdges[-(1:curSize), ]))]
      curEntries = cbind(rep(pos + (1:curNConst), each = curSize), allPermEdges)
      if (sparse) {
        Mat[cnt + (1:curNEntry), ] = curEntries
      } else {
        Mat[curEntries] = 1
      }
      cnt = cnt + curNEntry
      pos = pos + curNConst
    }
  }
  if (sparse) {
    Mat = sparseMatrix(i = Mat[, 1], j = Mat[, 2], x = 1, dims = c(numConstraints, numVars), index1 = TRUE)
  }
  output = list(mat = Mat, dir = Dir, rhs = Rhs, map = Map)
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

### This function identifies, among a list of graphs together with their bounds, a set of non-subsumed ones
### We say that a bound is subsumed if it is for a subgraph of a larger graph and is implied by its bound.
### For instance, the bound X(K_4 - e) >= 2 is subsumed by the bound X(K_4) >= 3 because of three criteria:
### a) K_4 - e is a subgraph of K_4 b) it has 1 less edge than K_4 and c) its bound is 1 less than for K_4.
### Similarly, the bound X(K_3 + e) >= 1 is subsumed by the bound X(K_3) >= 1 because of three criteria:
### a) K_3 + e is a supergraph of K_3 b) it has 1 more edge than K_3 and c) the bound is the same for K_3.
### The opposite reasoning applies to the upper bounds; note that only comparable graphs can be subsuming.
### Note that, for the purposes of an (I)LP a subsumed bound can only be tight if the subsuming one is too.
identifySubsumedBounds = function(allGraphs, boundTab, graphInfo, lower = TRUE) {
  boundTab %<>%
    filter(direction == ifelse(lower, "G", "L")) %>%
    group_by(graph) %>%
    arrange(bound * ifelse(lower, -1, 1)) %>%
    mutate(pos = row_number()) %>%
    mutate(subsumed = subsumed | (pos != 1)) %>%
    ungroup() %>%
    select(-pos)
  preSubsumed = boundTab %>%
    filter(subsumed)
  boundTab %<>%
    filter(!subsumed)
  L = nrow(boundTab)
  if (L == 1) { 
    output = bind_rows(preSubsumed, boundTab)
    return(output)
  }
  fullTab = boundTab %>%
    inner_join(graphInfo, by = join_by(graph == index)) %>%
    arrange(order, size) %>%
    mutate_at("bound", as.integer)
  print(paste("There are", L, ifelse(lower, "lower", "upper"), "bounds to process"))
  for (index in 1:L) {
    if (index %% 100 == 0) { print(index) }
    curRow      = fullTab %>% slice(index)
    curGraph    = allGraphs[[curRow$graph]]
    if (lower) {
      subsumedSuperCandidates = fullTab %>%
        filter(number != curRow$number & order >= curRow$order & size >= curRow$size & (bound - curRow$bound >= size - curRow$size))
      subsumedSubCandidates   = fullTab %>%
        filter(number != curRow$number & order <= curRow$order & size <= curRow$size & (bound - curRow$bound >= 0))
    } else {
      subsumedSuperCandidates = fullTab %>%
        filter(number != curRow$number & order >= curRow$order & size >= curRow$size & (curRow$bound - bound >= 0))
      subsumedSubCandidates   = fullTab %>%
        filter(number != curRow$number & order <= curRow$order & size <= curRow$size & (curRow$bound - bound >= curRow$size - size))
    }
    subsumedBySuper = sapply(allGraphs[subsumedSuperCandidates$graph], function(x) { subgraph_isomorphic(curGraph, x) })
    subsumedBySub   = sapply(allGraphs[  subsumedSubCandidates$graph], function(x) { subgraph_isomorphic(x, curGraph) })
    fullTab$subsumed[index] = (any(subsumedBySuper) || any(subsumedBySub))
  }
  fullTab %<>% select(number, support, graph, direction, bound, round, subsumed)
  output = bind_rows(preSubsumed, fullTab)
  output
}

### This function constructs the proof of a lower bound over the Ramsey polytope
### It makes use of the original LP as well as its dual variables at optimality
### If short = TRUE, only constructs a short proof version (with bound pointers)
constructProof = function(LP, dualVariables, lower = TRUE, short = FALSE) {
  goodRows   = which(!near(dualVariables, 0))
  dualValues = dualVariables[goodRows]
  dualValues = dualValues/min(abs(dualValues))
  dir        = LP$dir[goodRows]
  if (short) {
    relInfo   = LP$map[goodRows, , drop = FALSE]
    proof     = tibble(coeff = dualValues, constr = goodRows, dir = dir, bound = relInfo[,1], version = relInfo[,2], combo = relInfo[,3])
  } else {
    subMatrix = as.matrix(LP$mat[goodRows, , drop = FALSE])
    rhs       = LP$rhs[goodRows]
    proof     = cbind(coeff = dualValues, subMatrix, rhs = rhs)
    flipRows  = which(dir == ifelse(lower, "L", "G"))
    proof[flipRows, ] %<>% magrittr::multiply_by(-1)
  }
  proof
}

### This function computes the list of unique graphs isomorphic to a given graph
### If pairsOnly = FALSE, returns a 2m x p matrix (m: graph size, p: orbit size)
### If pairsOnly = TRUE,  returns a list of nC2 adjacency matrices, one per pair
permuteGraph = function(Graph, pairsOnly = FALSE) {
  if (pairsOnly) {
    Graph %<>% get.adjacency(type = "both", sparse = FALSE)
    n = nrow(Graph)
    stopifnot(ncol(Graph) == n)
    allPairs = combn2(1:n)
    allPerms = lapply(1:nrow(allPairs), function(x) { cur = allPairs[x,]; Graph[rev(cur),] = Graph[cur,]; Graph[, rev(cur)] = Graph[, cur]; Graph })
  } else {
    n = gorder(Graph)
    autoGens = automorphism_group(Graph)
    L = length(autoGens)
    autoSize = as.integer(automorphisms(Graph)$group_size)
    autoComp = sapply(autoGens, function(x) { min(which(x != 1:n)) })
    autoNext = sapply(1:L, function(x) { autoGens[[x]][autoComp[x]] })
    allPerms = matrix(unlist(permn(n)), nrow = n)
    for (ind in 1:L) {
      allPerms = allPerms[, allPerms[autoComp[ind], ] < allPerms[autoNext[ind], ], drop = FALSE]
    }
    numPerms = ncol(allPerms)
    expectedNumPerms = factorial(n)/autoSize
    Edges = as_edgelist(Graph)
    tailEdges = allPerms[Edges[,1], , drop = FALSE]
    headEdges = allPerms[Edges[,2], , drop = FALSE]
    swapPos   = which(tailEdges > headEdges)
    tempEdges = tailEdges[swapPos]
    tailEdges[swapPos] = headEdges[swapPos]
    headEdges[swapPos] = tempEdges
    if (numPerms == expectedNumPerms) { ### There are no duplicates, continue
      allPerms = rbind(tailEdges, headEdges)
    } else { ### There are duplicates; radix-sort adjacency lists and eliminate
      m = gsize(Graph)
      firstOrder  = cbind(as.vector(apply(headEdges, 2, order)), rep(1:numPerms, each = m))
      tailEdges   = matrix(tailEdges[firstOrder], ncol = numPerms)
      headEdges   = matrix(headEdges[firstOrder], ncol = numPerms)
      secondOrder = cbind(as.vector(apply(tailEdges, 2, order)), rep(1:numPerms, each = m))
      tailEdges   = matrix(tailEdges[secondOrder], ncol = numPerms)
      headEdges   = matrix(headEdges[secondOrder], ncol = numPerms)
      allPerms    = rbind(tailEdges, headEdges)
      allPerms    = allPerms[, !duplicated(t(allPerms)), drop = FALSE]
      stopifnot(ncol(allPerms) == expectedNumPerms)
    }
  }
  allPerms
}

### This function creates a matrix M with M[i,j] being the edge position of (ij)
### The edges are ordered as follows: 12, 13, ..., 1N, 23, ..., 2N, ..., (N-1)N
createPositionMatrix = function(numVertices) {
  numPos = choose(numVertices, 2)
  allPos = combn2(1:numVertices)
  outMat = matrix(0, numVertices, numVertices)
  outMat[allPos] = 1:numPos
  outMat
}

### This function traces back through a list of short proofs to obtain contradictions
reconstructContradictions = function(allProofs, boundTab) {
  lastIter = max(boundTab$iter)
  lastProofs = boundTab %>%
    filter(iter == lastIter) %>%
    group_by(graph) %>%
    mutate(N = n()) %>%
    filter(N == 2) %>%
    mutate(lower = max(bound * (direction == "G")), upper = max(bound * (direction == "L"))) %>%
    filter(upper < lower) %>%
    select(-N, -lower, -upper) %>%
    ungroup
  numContradictions = nrow(lastProofs)/2
  allContradictions = vector("list", numContradictions)
  for (index in 1:numContradictions) {
    curLastProof = lastProofs %>%
      slice(2 * (index - 1) + (1:2))
    curIndices = curLastProof$number
    allIndices = curIndices
    contradiction = vector("list", lastIter)
    pos = 1
    contradiction[[pos]] = curLastProof
    while (any(curIndices > 2)) {
      pos = pos + 1
      curIndices = setdiff(curIndices, 1:2)
      if (length(curIndices) > 0) {
        curProofs  = allProofs[curIndices]
        curIndices = lapply(curProofs, function(x) {x$bound}) %>%
          unlist() %>%
          c() %>%
          unique()
        contradiction[[pos]] = boundTab %>%
          slice(curIndices)
        allIndices %<>% c(curIndices)
      }
    }
    contradiction %<>% rev
    contradiction = contradiction[!sapply(contradiction, is.null)]
    contradiction %<>% 
      do.call(bind_rows, .) %>%
      distinct() %>%
      arrange(iter) %>%
      group_by(iter) %>%
      mutate(step = cur_group_id()) %>%
      ungroup()
    allContradictions[[index]] = contradiction
  }
  allContradictions
}

prepareBestProof = function(contradictions, listOfGraphs, labels = LETTERS, specialLabels = c("Z", "X"), fname = "ShortProofR3S4.RData") {
  numGraphs = sapply(contradictions, function(x) { n_distinct(x$graph) })
  numBounds = sapply(contradictions, nrow)
  numSteps  = sapply(contradictions, function(x) { n_distinct(x$iter) })
  contraTab = tibble(index = 1:length(contradictions), use = numGraphs, bounds = numBounds, steps = numSteps) %>%
    arrange(use, bounds, steps)
  bestProof = contraTab$index[1]
  bestContradiction = contradictions[[bestProof]]
  bestSize = numGraphs[bestProof]
  stopifnot(length(specialLabels) == 2)
  usableLabels = setdiff(labels, specialLabels)
  if (bestSize > length(usableLabels) + 2) { print("Not enough labels!"); return() }
  uGraphs = sort(unique(bestContradiction$graph))
  firstGraphs = bestContradiction %>%
    slice(1:2) %>%
    pull(graph)
  nameTab = tibble(index = setdiff(uGraphs, firstGraphs), label = usableLabels[1:(bestSize - 2)]) %>%
    bind_rows(tibble(index = firstGraphs, label = specialLabels)) %>%
    arrange(index)
  bestContradiction %<>%
    inner_join(nameTab, by = c("graph" = "index"))
  miniContradiction = bestContradiction %>%
    filter(step > 1) %>%
    mutate(summary = paste0("X[", label, "] \\", tolower(direction), "eq ", bound)) %>%
    arrange(step, label) %>%
    group_by(step) %>%
    mutate(fullSummary = paste0("$", paste(summary, collapse = ", "), "$")) %>%
    slice(1) %>%
    ungroup()
  allGraphObj = sapply(listOfGraphs, function(x) {graph_from_edgelist(x, directed = FALSE)})
  for (ind in 1:nrow(nameTab)) {
    curRow = nameTab %>%
      slice(ind)
    curGraph = allGraphObj[[curRow$index]]
    curLabel = curRow$label
    pdf(paste0(curLabel, ".pdf"))
    plot(curGraph)
    dev.off()
  }
  save(bestContradiction, file = fname)
  bestContradiction
}