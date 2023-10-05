library(combinat)
library(igraph)
library(cplexAPI)
library(magrittr)
library(Matrix)
library(slam)
library(tidyverse)

CPLEX_DIR        = "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/"

setwd("/Users/lchindelevitch/Downloads/NonPriority/Conjectures/RamseyNumbers")

### TODO: Consider implementing a version that looks for cyclic orientations too; the idea is to only constrain delta vectors which admit a bipartition!

### This function extracts connected graphs from the files provided by B. McKay
### Make sure that the conversion is done with the -el0o1 option in showg_mac64
### If extremeOnly = TRUE, only keeps the graphs with sizes within 2 of min.
getGraphs = function(numVerts = 2, extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE) {
  if (!treesOnly && !twoTreesOnly) {
    fn = paste0("Graphs", ifelse(extremeOnly, "Extreme", ""), numVerts, ".RData")
  } else {
    if (treesOnly) {
      fn = paste0("Trees", numVerts, ".RData")
    } else {
      fn = paste0("TwoTrees", numVerts, ".RData")
    }
  }
  if (!file.exists(fn)) {
    curText = readLines(paste0(ifelse(treesOnly, "trees", ifelse(twoTreesOnly, "TwoTrees", "graph")), numVerts, ifelse(extremeOnly, "x", "c"), ".txt"))
    L = length(curText)
    goodLines = 4 * (1:(L/4))
    numGraphs = L/4
    numEdges = curText[goodLines - 1] %>%
      str_split_fixed(" ", n = 2) %>%
      magrittr::extract(, 2) %>%
      as.integer
    curText = curText[goodLines]
    curText = split(curText, numEdges)
    numByEdge = sapply(curText, length)
    numCounts = length(curText)
    allG = vector("list", numCounts)
    pos = 1
    for (index in 1:numCounts) {
      curNumber = numByEdge[index]
      curCounts = as.integer(names(curText)[index])
      print(paste("Processing the", curNumber, "graphs with", curCounts, "edges"))
      if (!extremeOnly || (curCounts <= numVerts) || (curCounts >= choose(numVerts, 2) - 1)) {
        curX = curText[[index]] %>%
          str_split("  ") %>%
          unlist() %>%
          str_split_fixed(" ", n = 2)
        mode(curX) = "integer"
        tails = curX[,1] %>%
          matrix(nrow = curCounts)
        heads = curX[,2] %>%
          matrix(nrow = curCounts)
        curX = rbind(tails, heads)
        allG[[pos]] = as.list(as.data.frame(curX)) %>%
          set_names(NULL)
        pos = pos + 1
      }
    }
    allG = allG[1:(pos - 1)]
    allG = do.call(c, allG)
    save(allG, file = fn)
  }
  else {
    e = new.env()
    load(fn, envir = e)
    allG = get("allG", envir = e)
  }
  allG
}

findMinimalGraphSet = function(n = 9L, r = 3L, s = 4L, maxOrder = 6L, symRS = FALSE, shortProofs = TRUE, 
                               minAuto = NULL, extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, inds = NULL) {
  lastSolution = NULL
  inds = sort(setdiff(inds, (1:ifelse(r == s, 1, 2))))
  L = length(inds)
  print(paste("There are", L, "indices to process"))
  pos = L
  while (pos >= 1) {
    print(pos)
    altInds = inds[-pos]
    altRes  = iterateLPProof(n = n, r = r, s = s, maxOrder = maxOrder, symRS = symRS, shortProofs = shortProofs,
                             minAuto = minAuto, extremeOnly = extremeOnly, treesOnly = treesOnly, twoTreesOnly = twoTreesOnly, inds = altInds)
    if (length(altRes) > 1) {
      print(paste("Removing", inds[pos], "leaves the contradiction valid"))
      lastSolution = altRes
      bestSolution = prepareBestProof(altRes[[4]], altRes[[3]], labels = c(LETTERS, letters), plotGraphs = FALSE)
      inds = altInds[sort(unique(setdiff(bestSolution$graph, (1:ifelse(r == s, 1, 2)))))]
      pos = length(inds)
    } else {
      pos = pos - 1
    }
  }
  output = list(inds = inds, last = lastSolution)
  output
}

### This function tries to iteratively construct a proof that Ramsey(r,s) <= n.
### Search parameters (apply to the LPs and the corresponding proofs):
### maxOrder determines the largest order of the graphs that will be considered.
### If symRS = TRUE, only lower bounds are computed, the upper bounds being symmetric.
### If shortProofs = TRUE, only produces the shortened (not full) version of proofs.
### Graph selection parameters (Kr and Ks are considered essential and put 1st/2nd):
### If minAuto is not NULL, it is a lower bound on the automorphism group size.
### If extremeOnly = TRUE, only keeps the graphs with sizes within 1 of min/max.
### If treesOnly = TRUE, only keeps the graphs with size equal to the order - 1.
### If twoTreesOnly = TRUE, only keeps the graphs with size equal to the 2 * order - 3.
### If inds is not NULL, only keeps the non-essential graphs whose numbers are in inds.
iterateLPProof = function(n = 6L, r = 3L, s = r, maxOrder = 2 * n / 3, symRS = (r == s), shortProofs = TRUE,
                          minAuto = NULL, extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, inds = NULL) {
  allGraphs = list(as.vector(combn2(1:r)))
  if (r !=s ) {
    allGraphs %<>% c(list(as.vector(combn2(1:s))))
  }
  for (gsize in 2:maxOrder) {
    allGraphs %<>% c(getGraphs(gsize, extremeOnly = extremeOnly, treesOnly = treesOnly, twoTreesOnly = twoTreesOnly))
  }
  numAutos    = sapply(allGraphs, function(x) { as.integer(automorphisms(graph_from_edgelist(matrix(x, ncol = 2), directed = FALSE))$group_size) })
  graphTab    = tibble(size = sapply(allGraphs, length)/2, order = sapply(allGraphs, max), n_auto = numAutos, orbit = factorial(order)/n_auto)
  repeatGraphs = graphTab %>%
    rowid_to_column(var = "index") %>%
    filter(order %in% c(r, s) & size == choose(order, 2) & index > ifelse(r == s, 1, 2)) %>%
    pull(index)
  if (length(repeatGraphs) > 0) {
    allGraphs %<>% magrittr::extract(-repeatGraphs)
    graphTab  %<>% slice(-repeatGraphs)
  }
  goodGraphs   = (graphTab$order <= maxOrder)
  if (extremeOnly) { 
    goodGraphs = (goodGraphs & (graphTab$size <= graphTab$order) | graphTab$size >= (choose(graphTab$order, 2) - 1)) 
  }
  goodGraphs   = (goodGraphs & (graphTab$n_auto >= ifelse(is.null(minAuto), 1, minAuto))) 
  goodGraphs   = sort(unique(c(1:ifelse(r == s, 1, 2), which(goodGraphs))))
  if (!is.null(inds)) {
    inds = sort(unique(c(1:ifelse(r == s, 1, 2), inds)))
    goodGraphs %<>% magrittr::extract(inds)
  }
  allGraphs %<>% magrittr::extract(goodGraphs)
  graphTab  %<>% slice(goodGraphs) %>%
    rowid_to_column(var = "index")
  numGraphs = length(goodGraphs)
  allGraphAdj = vector("list", numGraphs)
  boundTab = tibble(number = 1:2, support = c(r, s), graph = c(1, ifelse(r == s, 1, 2)), direction = c("G", "L"), bound = c(1, choose(s, 2) - 1), round = 0L, subsumed = FALSE, size = choose(c(r, s), 2))
  boundIndex = 3
  proofs     = vector("list", n^2)
  curSupport = min(r, s) + 1
  curOrder   = curSupport
  curRound   = 0
  while (curOrder <= maxOrder) {
    curGraphTab = graphTab %>% 
      filter(order >= min(r, s) & order <= curOrder)
    print(paste("Exploring graphs on up to", curOrder, "vertices with support size", curSupport))
    lowerBounds = identifySubsumedBounds(allGraphs, boundTab, curGraphTab, lower = TRUE)
    upperBounds = identifySubsumedBounds(allGraphs, boundTab, curGraphTab, lower = FALSE)
    boundTab    = bind_rows(lowerBounds, upperBounds)
    curRound    = curRound + 1
    finalResult  = findNextGraph(numVertices = curSupport, allGraphs = allGraphs, allGraphAdj = allGraphAdj, boundTab = boundTab, 
                                 graphInfo = curGraphTab, nRound = paste0("R", r, "S", s, "I", curRound), symRS = symRS, shortProofs = shortProofs)
    fullResult  = finalResult[[1]]
    allGraphAdj = finalResult[[2]]
    contradict = FALSE
    goodPieces = !sapply(fullResult, is.null)
    improved = any(goodPieces)
    if (improved) {
      fullResult = fullResult[goodPieces]
      for (ind in 1:length(fullResult)) {
        curResult = fullResult[[ind]]
        curGraph  = curResult$graph
        curSize = curGraphTab %>% 
          filter(index == curGraph) %>% 
          pull(size)
        if (!is.null(curResult$proofL)) {
          proofs[[boundIndex]] = curResult$proofL
          boundTab %<>% 
            bind_rows(tibble(number = boundIndex, support = curSupport, graph = curGraph, direction = "G", bound = curResult$boundL, round = curRound, subsumed = FALSE, size = curSize))
          boundIndex = boundIndex + 1
        }
        if (!is.null(curResult$proofU)) {
          proofs[[boundIndex]] = curResult$proofU
          boundTab %<>% 
            bind_rows(tibble(number = boundIndex, support = curSupport, graph = curGraph, direction = "L", bound = curResult$boundU, round = curRound, subsumed = FALSE, size = curSize))
          boundIndex = boundIndex + 1
        }
        contradict = (contradict || (curResult$boundL > curResult$boundU) || (curResult$boundL == curSize) || (curResult$boundU == 0))
      }
    } else {
      curSupport = curSupport + 1
      if (curOrder < maxOrder) {
        curOrder = curSupport
      } else {
        if (curSupport > n) {
          print("Not enough evidence for a contradiction!")
          return(FALSE)
        }
      }
    }
    if (contradict) {
      print("Found a contradiction!")
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
  firstGraphs = lapply(allGraphs, function(x) { matrix(x, ncol = 2) })
  output = list(proofs = proofs, boundTab = boundTab, graphs = firstGraphs, contradictions = contradictions)
  output
}

### This function finds bounds on the number of red edges in a graph
### The Ramsey polytope is constructed over the specified order, numVertices; known bounds are provided via boundTab
### The output will contain information on all of the input graphs provided (allGraphs) that yield a non-trivial bound
### The allGraphAdj variable contains all the distinct permuted versions of each graph's adjacency list in the same order
### that also contain relevant information in graphInfo (in particular, whose indices are contained in its first column).
findNextGraph = function(numVertices, allGraphs, allGraphAdj, boundTab, graphInfo, nRound = 0, symRS = FALSE, shortProofs = TRUE) {
  L = nrow(graphInfo)
  prep = prepareRamseyLP(numVertices = numVertices, allGraphs = allGraphs, allGraphAdj = allGraphAdj, graphBounds = boundTab, graphInfo = graphInfo, sparse = TRUE)
  fname    = paste0("TempFile", nRound, ".lp")
  numVars  = ncol(prep$mat)
  numConst = length(prep$dir)
  numNonZeros = tail(prep$mat@p,1)
  CN = paste0("X", 1:numVars)
  ### Writing the model (with no objective function!) into an LP file with cplexAPI
  envir = cplexAPI::openEnvCPLEX()
  model = cplexAPI::initProbCPLEX(envir, pname = "Ramsey Polytope Optimisation")
  cplexAPI::copyLpwNamesCPLEX(env = envir, lp = model, nCols = numVars, nRows = numConst, lpdir = CPX_MIN, objf = rep(0, numVars), rhs = prep$rhs, sense = prep$dir, 
  matbeg = prep$mat@p[1:numVars], matcnt = diff(prep$mat@p), matind = prep$mat@i, matval = rep(1, numNonZeros), lb = rep(0, numVars), ub = rep(1, numVars), cnames = CN)
  cplexAPI::writeProbCPLEX(envir, model, fname = fname)
  print(paste("There are", numConst, "constraints over", numVars, "variables"))
  output      = vector("list", L)
  auxMatrix   = createPositionMatrix(numVertices)
  print(paste("There are", L, "graphs to process"))
  initObjective = rep(0, numVars)
  for (ind in 1:L) {
    print(ind)
    # curOutFiles = outputFiles[ind * 2 - (1:0)]
    curRow     = graphInfo %>% slice(ind)
    curIndex   = curRow$index
    curGraph   = allGraphs[[curIndex]]
    curSize    = length(curGraph)/2
    curEdges   = auxMatrix[matrix(curGraph, ncol = 2)]
    curObjective = initObjective
    curObjective[curEdges] = 1
    cplexAPI::setIntParmCPLEX(envir, cplexAPI::CPXPARAM_Preprocessing_Presolve, 0)
    cplexAPI::chgObjCPLEX(envir, model, ncols = numVars, ind = 0:(numVars - 1), val = curObjective)
    cplexAPI::dualoptCPLEX(envir, model)
    curSolution = cplexAPI::solutionCPLEX(envir, model)
    curObjL = curSolution$objval
    boundL  = as.integer(ifelse(near(curObjL, round(curObjL)), round(curObjL), ceiling(curObjL)))
    proofL  = NULL
    prevLowerBound = boundTab %>% 
      filter(graph == curIndex & direction == "G") %>% 
      pull(bound)
    if (boundL > max(c(0, prevLowerBound))) {
      curDualsL  = curSolution$pi
      if (!(all(near(curDualsL[1:numConst], 0)))) {
        proofL   = constructProof(prep, curDualsL[1:numConst], lower = TRUE, shortProofs = shortProofs)
      }
    }
    if (symRS) {
      boundU = curSize - boundL
      proofU = proofL
    } else {
      cplexAPI::setObjDirCPLEX(envir, model, CPX_MAX)
      cplexAPI::dualoptCPLEX(envir, model)
      altSolution = cplexAPI::solutionCPLEX(envir, model)
      curObjU = altSolution$objval
      boundU  = as.integer(ifelse(near(curObjU, round(curObjU)), round(curObjU), floor(curObjU)))
      proofU  = NULL
      prevUpperBound = boundTab %>% 
        filter(graph == curIndex & direction == "L") %>% 
        pull(bound)
      if (boundU < min(c(curSize, prevUpperBound))) {
        curDualsU  = altSolution$pi
        if (!(all(near(curDualsU[1:numConst], 0)))) {
          proofU   = constructProof(prep, curDualsU[1:numConst], lower = FALSE, shortProofs = shortProofs)
        }
      }
      cplexAPI::setObjDirCPLEX(envir, model, CPX_MIN)
    }
    curResult = NULL
    if (!is.null(proofL) || !is.null(proofU)) {
      curResult = list(graph = curIndex, proofL = proofL, proofU = proofU, boundL = boundL, boundU = boundU)
    }
    output[[ind]] = curResult
  }
  cplexAPI::delProbCPLEX(envir, model)
  cplexAPI::closeEnvCPLEX(envir)
  finalOutput = list(output = output, adj = prep$adj)
  finalOutput
}

### This function prepares to optimise over the Ramsey polytope with specified bounds
prepareRamseyLP = function(numVertices, allGraphs, allGraphAdj, graphBounds, graphInfo, sparse = TRUE) {
  numVars = choose(numVertices, 2)
  auxMatrix = createPositionMatrix(numVertices = numVertices)
  graphBounds %<>%
    filter(!subsumed) %>%
    inner_join(graphInfo, by = join_by(graph == index, size == size)) %>%
    arrange(graph, direction) %>%
    mutate(numOpts = choose(numVertices, order), numCombinations = orbit * numOpts, fullCombinations = numCombinations * size)
  numConstraints = sum(graphBounds$numCombinations)
  Dir = rep("", numConstraints)
  Rhs = rep(0,  numConstraints)
  Szs = rep(0, numConstraints)
  Map = matrix(0,  numConstraints, 3)
  Mat = matrix(0, ifelse(sparse, sum(graphBounds$fullCombinations), numConstraints), ifelse(sparse, 2, numVars))
  pos = 0
  cnt = 0
  uGraphs = unique(graphBounds$graph)
  print(paste("There are", length(uGraphs), "graphs to prepare for the LP"))
  for (iter in 1:length(uGraphs)) {
    if (iter %% 10 == 0) { print(iter) }
    ind = uGraphs[iter]
    if (is.null(allGraphAdj[[ind]])) {
      allGraphAdj[[ind]] = permuteGraph(allGraphs[[ind]])
      gc()
    }
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
    Szs[pos + (1:curNFull)]   = rep(curSize,              curCombos)
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
  output = list(mat = Mat, dir = Dir, rhs = Rhs, szs = Szs, map = Map, adj = allGraphAdj)
  output
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
    inner_join(graphInfo, by = join_by(graph == index, size == size)) %>%
    arrange(order, size) %>%
    mutate_at("bound", as.integer)
  print(paste("There are", L, ifelse(lower, "lower", "upper"), "bounds to process"))
  for (index in 1:L) {
    if (index %% 100 == 0) { print(index) }
    curRow      = fullTab %>% slice(index)
    curGraph    = graph_from_edgelist(matrix(allGraphs[[curRow$graph]], ncol = 2), directed = FALSE)
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
    subsumedBySuper = sapply(allGraphs[subsumedSuperCandidates$graph], function(x) { subgraph_isomorphic(curGraph, graph_from_edgelist(matrix(x, ncol = 2), directed = FALSE)) })
    subsumedBySub   = sapply(allGraphs[  subsumedSubCandidates$graph], function(x) { subgraph_isomorphic(graph_from_edgelist(matrix(x, ncol = 2), directed = FALSE), curGraph) })
    fullTab$subsumed[index] = (any(subsumedBySuper) || any(subsumedBySub))
  }
  fullTab %<>% select(number, support, graph, direction, bound, round, subsumed, size)
  output = bind_rows(preSubsumed, fullTab)
  output
}

### This function constructs the proof of a lower bound over the Ramsey polytope
### It makes use of the original LP as well as its dual variables at optimality
### If lower = FALSE the funciton constructs the proof of an upper bound instead
### If shortProofs = TRUE, only constructs a short proof version (with bound pointers)
constructProof = function(LP, dualVariables, lower = TRUE, shortProofs = FALSE) {
  goodRows   = which(!near(dualVariables, 0))
  dualValues = dualVariables[goodRows]
  dualValues = dualValues/min(abs(dualValues))
  dir        = LP$dir[goodRows]
  if (shortProofs) {
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
### Returns a 2m x p matrix (m: graph size, p: orbit size)
permuteGraph = function(Graph) {
  Graph  = matrix(Graph, ncol = 2)
  n = max(Graph)
  iGraph = graph_from_edgelist(Graph, directed = FALSE)
  autoGens = automorphism_group(iGraph)
  L = length(autoGens)
  autoSize = as.integer(automorphisms(iGraph)$group_size)
  allPerms = permn(n)
  if (L > 0) {
    autoComp = sapply(autoGens, function(x) { min(which(x != 1:n)) })
    autoNext = sapply(1:L, function(x) { autoGens[[x]][autoComp[x]] })
    for (ind in 1:L) {
      allPerms = allPerms[map_int(allPerms, ~{nth(., autoComp[ind])}) < map_int(allPerms, ~{nth(., autoNext[ind])})]
    }
  }
  allPerms = matrix(unlist(allPerms), nrow = n)
  numPerms = ncol(allPerms)
  expectedNumPerms = factorial(n)/autoSize
  tailEdges = allPerms[Graph[,1], , drop = FALSE]
  headEdges = allPerms[Graph[,2], , drop = FALSE]
  swapPos   = which(tailEdges > headEdges)
  tempEdges = tailEdges[swapPos]
  tailEdges[swapPos] = headEdges[swapPos]
  headEdges[swapPos] = tempEdges
  if (numPerms == expectedNumPerms) { ### There are no duplicates, continue
    allPerms = rbind(tailEdges, headEdges)
  } else { ### There are duplicates; radix-sort adjacency lists and eliminate
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
    filter(N == 2)
  if (nrow(lastProofs) > 0) {
    lastProofs %<>%
      mutate(lower = max(bound * (direction == "G")), upper = max(bound * (direction == "L"))) %>%
      filter(upper < lower) %>%
      select(-N, -lower, -upper) %>%
      ungroup
  } else {
    lastProofs %<>% 
      select(-N) %>%
      ungroup()
  }
  extraProofs = boundTab %>%
    filter(iter == lastIter) %>%
    filter((direction == "G" & bound == size) | (direction == "L" & bound == 0))
  lastProofs %<>%
    bind_rows(extraProofs)
  uGraphs = sort(unique(lastProofs$graph))
  numContradictions = length(uGraphs)
  allContradictions = vector("list", numContradictions)
  for (index in 1:numContradictions) {
    curGraph = uGraphs[index]
    curLastProof = lastProofs %>%
      filter(graph == curGraph)
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

### This function prepares the proof of a contradiction using the fewest graphs and iterations.
### The first input is the list of contradictions, the second, the list of all the graphs used.
### The labels are used to create a file with each relevant graph; the special labels are used
### for the two starting graphs (which can also be identical); the best proof is saved in fname.
prepareBestProof = function(contradictions, listOfGraphs, labels = LETTERS, specialLabels = c("Z", "X"), fname = "ShortProofR3S4.RData", plotGraphs = FALSE) {
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
    inner_join(nameTab, by = join_by(graph == index))
  miniContradiction = bestContradiction %>%
    filter(step > 1) %>%
    mutate(summary = paste0("X[", label, "] \\", tolower(direction), "eq ", bound)) %>%
    arrange(step, label) %>%
    group_by(step) %>%
    mutate(fullSummary = paste0("$", paste(summary, collapse = ", "), "$")) %>%
    slice(1) %>%
    ungroup()
  if (plotGraphs) {
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
  }
  save(bestContradiction, miniContradiction, file = fname)
  bestContradiction
}

### This function parses a written proof in LaTeX format, providing the totals. 
### Note: replace all \ with \\ in the LaTeX source before passing it as input!
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
