gc()

library(combinat)
library(igraph)
library(cplexAPI)
library(magrittr)
library(Matrix)
library(slam)
library(tidyverse)

CPLEX_DIR        = "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/"

MAX_SIZE = 10L
print("Initializing permutations")
for (ind in 1:MAX_SIZE) {
  print(ind)
  assign(paste0("PERMS", ind), permn(ind), envir = .GlobalEnv)
}

### GOAL: Check whether the constraints always form a set of disjoint trees, and double-check with Knuth's formula - in fact, it may be used to generate the embeddings!

setwd("/Users/lchindelevitch/Downloads/NonPriority/Conjectures/RamseyNumbers")

### TODO: Consider implementing a version that looks for cyclic orientations too; the idea is to only constrain delta vectors which admit a bipartition!

### This function extracts connected graphs from the files made by McKay's nauty
### Make sure that the conversion is done with the -el0o1 option in nauty::showg
### If extremeOnly = TRUE, only keeps the graphs with sizes within 2 of min.
getGraphs = function(numVerts = 2, extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, ETOnly = FALSE, partiteOnly = FALSE) {
  if (!treesOnly && !twoTreesOnly && !ETOnly && !partiteOnly) {
    fn = paste0("Graphs", ifelse(extremeOnly, "Extreme", ""), numVerts, ".RData")
  } else {
    if (treesOnly) {
      fn = paste0("Trees", numVerts, ".RData")
    }
    if (twoTreesOnly) {
      fn = paste0("TwoTrees", numVerts, ".RData")
    }
    if (ETOnly) {
      fn = paste0("ET", numVerts, "c.RData")
    }
    if (partiteOnly) {
      fn = paste0("Bipartite", numVerts, ".RData")
    }
  }
  if (!file.exists(fn)) {
    curText = readLines(paste0(ifelse(treesOnly, "trees", ifelse(twoTreesOnly, "TwoTrees", ifelse(ETOnly, "ET", ifelse(partiteOnly, "Bipartite", "graph")))), numVerts, ifelse(extremeOnly, "x", "c"), ".txt"))
    L = length(curText)
    goodLines = 4 * (1:(L/4))
    numGraphs = L/4
    numEdges = curText[goodLines - 1] %>%
      str_split_fixed(" ", n = 2) %>%
      magrittr::extract(, 2) %>%
      as.integer()
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

findMinimalGraphSet = function(n = 9L, r = 3L, s = 4L, maxOrder = 6L, symRS = FALSE, shortProofs = TRUE, lowerOnly = FALSE, minAuto = NULL, eps = NA,
                               extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, ETOnly = FALSE, partiteOnly = FALSE, inds = NULL) {
  lastSolution = NULL
  endpoint = ifelse(r == s, 1, 2)
  inds = sort(setdiff(inds, (1:endpoint)))
  L = length(inds)
  print(paste("There are", L, "indices to process"))
  pos = L
  while (pos >= 1) {
    print(pos)
    altInds = inds[-pos]
    altRes  = iterateLPProof(n = n, r = r, s = s, maxOrder = maxOrder, symRS = symRS, shortProofs = shortProofs, lowerOnly = lowerOnly, minAuto = minAuto, eps = eps,
                             extremeOnly = extremeOnly, treesOnly = treesOnly, twoTreesOnly = twoTreesOnly, ETOnly = ETOnly, partiteOnly = partiteOnly, inds = altInds)
    if (length(altRes) > 1) {
      print(paste("Removing", inds[pos], "leaves the contradiction valid"))
      lastSolution = altRes
      bestSolution = prepareBestProof(altRes[[4]], altRes[[3]], labels = c(LETTERS, letters), plotGraphs = FALSE)
      inds = altInds[sort(unique(setdiff(bestSolution$graph, (1:endpoint)))) - endpoint]
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
### If lowerOnly = TRUE, no upper, bounds are computed (this is regardless of symRS).
### Graph selection parameters (Kr and Ks are considered essential and put 1st/2nd):
### If minAuto is not NULL, it is a lower bound on the automorphism group size.
### If eps is not NA, the current bound needs to be tightened by at least eps to count.
### If extremeOnly = TRUE, only keeps the graphs with sizes within 1 of min/max.
### If treesOnly = TRUE, only keeps the graphs with size equal to the order - 1.
### If twoTreesOnly = TRUE, only keeps the 2-trees (size equal to 2 * order - 3).
### If ETOnly = TRUE, only keeps the edge-transitive graphs from the OEIS reference.
### If partiteOnly = TRUE, only keeps the bipartite graphs as generated by nauty.
### If completeOnly = TRUE, only keeps the complete and complete bipartite graphs.
### If inds is not NULL, only keeps the non-essential graphs whose numbers are in inds.
iterateLPProof = function(n = 6L, r = 3L, s = r, maxOrder = round(2 * n / 3), symRS = (r == s), shortProofs = TRUE, lowerOnly = FALSE, minAuto = NULL, eps = NA,
                          extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, ETOnly = FALSE, partiteOnly = FALSE, completeOnly = FALSE, inds = NULL) {
  allGraphs = list(as.vector(combn2(1:r)))
  if (r !=s ) {
    allGraphs %<>% c(list(as.vector(combn2(1:s))))
  }
  endpoint = ifelse(r == s, 1, 2)
  if (!completeOnly) {
    for (gsize in 2:maxOrder) {
      allGraphs %<>% c(getGraphs(gsize, extremeOnly = extremeOnly, treesOnly = treesOnly, twoTreesOnly = twoTreesOnly, ETOnly = ETOnly, partiteOnly = partiteOnly))
    }
    numAutos    = sapply(allGraphs, function(x) { as.integer(automorphisms(graph_from_edgelist(matrix(x, ncol = 2), directed = FALSE))$group_size) })
    graphTab    = tibble(size = sapply(allGraphs, length)/2, order = sapply(allGraphs, max), n_auto = numAutos, n_embed = factorial(order)/n_auto)
    repeatGraphs = graphTab %>%
      rowid_to_column(var = "index") %>%
      filter(order %in% c(r, s) & size == choose(order, 2) & index > endpoint) %>%
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
    goodGraphs   = sort(unique(c(1:endpoint, which(goodGraphs))))
  } else {
    numAutos = unique(c(factorial(r), factorial(s)))
    for (gsize in 2:maxOrder) {
      if (gsize != r && gsize != s) {
        allGraphs %<>% c(list(as.vector(combn2(1:gsize))))
        numAutos  %<>% c(factorial(gsize))
      }
      for (left in ceiling(gsize/2):(gsize - 1)) {
        right  = gsize - left
        nextGraph = hcube(c(left, right))
        nextGraph[,2] %<>% add(left)
        allGraphs %<>% c(list(as.vector(nextGraph)))
        numAutos  %<>% c(factorial(left) * factorial(right) * ifelse(left == right, 2, 1))
      }
    }
    graphTab = tibble(size = sapply(allGraphs, length)/2, order = sapply(allGraphs, max), n_auto = numAutos, n_embed = factorial(order)/n_auto)
    goodGraphs = 1:nrow(graphTab)
  }
  if (!is.null(inds)) {
    inds = sort(unique(c(1:endpoint, inds)))
    goodGraphs %<>% magrittr::extract(inds)
  }
  allGraphs %<>% magrittr::extract(goodGraphs)
  graphTab  %<>% slice(goodGraphs) %>%
    rowid_to_column(var = "index")
  numGraphs = length(goodGraphs)
  allGraphAdj = vector("list", numGraphs)
  allGraphOrb = vector("list", numGraphs)
  boundTab = tibble(number = 1:2, support = c(r, s), graph = c(1, endpoint), direction = c("G", "L"), bound = c(1, choose(s, 2) - 1), round = 0L, subsumed = FALSE, size = choose(c(r, s), 2))
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
    finalResult  = tightenBounds(numVertices = curSupport, allGraphs = allGraphs, allGraphAdj = allGraphAdj, allGraphOrb = allGraphOrb, boundTab = boundTab, 
                   graphInfo = curGraphTab, nRound = paste0("R", r, "S", s, "I", curRound), symRS = symRS, shortProofs = shortProofs, eps = eps, lowerOnly = lowerOnly)
    fullResult  = finalResult[[1]]
    allGraphAdj = finalResult[[2]]
    allGraphOrb = finalResult[[3]]
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
        contradict = (contradict || (curResult$boundL == curSize) || (!is.null(curResult$boundU) && (curResult$boundU == 0 || (curResult$boundL > curResult$boundU))))
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

### This function looks for tighter bounds on the number of red edges in a graph using currently available constraints.
### The Ramsey polytope is constructed over the specified order, numVertices; known bounds are provided via boundTab
### The output will contain information on all of the input graphs provided (allGraphs) that yield a non-trivial bound
### The allGraphAdj variable contains all the distinct permuted versions of each graph's adjacency list in the same order
### The allGraphOrb variable contains all the information about the orbits of each graph's automorphism group
### that also contain relevant information in graphInfo (in particular, whose indices are contained in its first column).
### If eps is specified (ie, not NA), the function looks for an improvement over the best existing bound by at least eps.
tightenBounds = function(numVertices, allGraphs, allGraphAdj, allGraphOrb, boundTab, graphInfo, nRound = 0, symRS = FALSE, shortProofs = TRUE, eps = NA, lowerOnly = FALSE) {
  L = nrow(graphInfo)
  prep = prepareRamseyLP(numVertices = numVertices, allGraphs = allGraphs, allGraphAdj = allGraphAdj, allGraphOrb = allGraphOrb, graphBounds = boundTab, 
                         graphInfo = graphInfo, sparse = TRUE)
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
  cplexAPI::setIntParmCPLEX(envir, cplexAPI::CPXPARAM_Preprocessing_Presolve, 0)
  print(paste("There are", numConst, "constraints over", numVars, "variables"))
  output      = vector("list", L)
  auxMatrix   = createPositionMatrix(numVertices, symmetric = TRUE)
  print(paste("There are", L, "graphs to process"))
  initObjective = rep(0, numVars)
  for (ind in 1:L) {
    print(ind)
    curRow     = graphInfo %>% slice(ind)
    curIndex   = curRow$index
    curOrder   = curRow$order
    curGraph   = allGraphs[[curIndex]]
    curSize    = length(curGraph)/2
    curEdges   = auxMatrix[matrix(curGraph, ncol = 2)]
    prevLowerBound = max(c(0,       boundTab %>% filter(graph == curIndex & direction == "G") %>% pull(bound)))
    prevUpperBound = min(c(curSize, boundTab %>% filter(graph == curIndex & direction == "L") %>% pull(bound)))
    proofL  = NULL
    proofU  = NULL
    boundL  = NULL
    boundU  = NULL
    curResult = NULL
    curObjective = initObjective
    redCNames = cplexAPI::getColNameCPLEX(envir, model, 0, numVars - 1) %>% str_remove_all("X") %>% as.integer()
    curObjective[match(curEdges, redCNames)] = 1
    cplexAPI::chgObjCPLEX(envir, model, ncols = numVars, ind = 0:(numVars - 1), val = curObjective)
    numExtraConst = 0
    if (!is.na(eps)) {
      cplexAPI::setDblParmCPLEX(envir, cplexAPI::CPXPARAM_Simplex_Limits_LowerObj, prevLowerBound + eps)
    }
    if (is.null(allGraphOrb[[curIndex]])) {
      allGraphOrb[[curIndex]] = getOrbits(curGraph)
    }
    curGraphOrb = allGraphOrb[[curIndex]]
    extraOrbits = curGraphOrb$orbitsE
    if (curOrder < numVertices) {
      numExtraVerts = numVertices - curOrder
      extraVerts = (curOrder + 1):numVertices
      if (numExtraVerts > 1) {
        extraOutsideOrbits = as.vector(combn2(extraVerts))
        extraOrbits = c(extraOrbits, list(extraOutsideOrbits))
      }
      vOrbits = curGraphOrb$orbitsV
      extraInsideOutsideOrbits = lapply(vOrbits, function(x) { as.vector(as.matrix(expand.grid(x, extraVerts))) })
      extraOrbits = c(extraOrbits, extraInsideOutsideOrbits)
    }
    extraOrbits = lapply(extraOrbits, function(x) { auxMatrix[matrix(x, ncol = 2)]} )
    extraOrbits = extraOrbits[sapply(extraOrbits, length) > 1]
    if (length(extraOrbits) > 0) {
      numExtraConst = sum(sapply(extraOrbits, length)) - length(extraOrbits)
      extraIndices = lapply(extraOrbits, function(x) { y = rep(x, each = 2); y = y[-c(1, length(y))]; y })
      cplexAPI::addRowsCPLEX(envir, model, ncols = 0, nrows = numExtraConst, nnz = 2 * numExtraConst, sense = rep("E", numExtraConst), rhs = rep(0, numExtraConst),
                           matbeg = 2 * (0:(numExtraConst - 1)), matind = match(unlist(extraIndices), redCNames) - 1, matval = rep(c(1,-1), numExtraConst))
    }
    ### cplexAPI::writeProbCPLEX(envir, model, fname = fname)
    cplexAPI::dualoptCPLEX(envir, model)
    curSolution = cplexAPI::solutionCPLEX(envir, model)
    if (curSolution$lpstat == cplexAPI::CPX_STAT_OPTIMAL) {
      curObjL = curSolution$objval
      boundL  = as.integer(ifelse(near(curObjL, round(curObjL)), round(curObjL), ceiling(curObjL)))
      if (boundL > prevLowerBound) {
        curDualsL  = curSolution$pi
        if (!(all(near(curDualsL, 0)))) {
          proofL   = constructProof(prep, curDualsL[1:numConst], lower = TRUE, shortProofs = shortProofs)
        }
      }
    }
    if (!is.na(eps)) {
      cplexAPI::setDblParmCPLEX(envir, cplexAPI::CPXPARAM_Simplex_Limits_LowerObj, 0)
    }
    if (symRS) {
      boundU = curSize - boundL
      proofU = proofL
    } else {
      if (!lowerOnly) {
        cplexAPI::setObjDirCPLEX(envir, model, CPX_MAX)
        if (!is.na(eps)) {
          cplexAPI::setDblParmCPLEX(envir, cplexAPI::CPXPARAM_Simplex_Limits_UpperObj, prevUpperBound - eps)
        }
        cplexAPI::dualoptCPLEX(envir, model)
        altSolution = cplexAPI::solutionCPLEX(envir, model)
        if (curSolution$lpstat == cplexAPI::CPX_STAT_OPTIMAL) {
          curObjU = altSolution$objval
          boundU  = as.integer(ifelse(near(curObjU, round(curObjU)), round(curObjU), floor(curObjU)))
          if (boundU < prevUpperBound) {
            curDualsU  = altSolution$pi
            if (!(all(near(curDualsU, 0)))) {
              proofU   = constructProof(prep, curDualsU[1:numConst], lower = FALSE, shortProofs = shortProofs)
            }
          }
        }
        if (!is.na(eps)) {
          cplexAPI::setDblParmCPLEX(envir, cplexAPI::CPXPARAM_Simplex_Limits_UpperObj, curSize)
        }
        cplexAPI::setObjDirCPLEX(envir, model, CPX_MIN)
      }
    }
    if (numExtraConst > 0) {
      cplexAPI::delRowsCPLEX(envir, model, numConst, numConst + numExtraConst - 1)
    }
    if (!is.null(proofL) || !is.null(proofU)) {
      curResult = list(graph = curIndex, proofL = proofL, proofU = proofU, boundL = boundL, boundU = boundU)
    }
    output[[ind]] = curResult
  }
  cplexAPI::delProbCPLEX(envir, model)
  cplexAPI::closeEnvCPLEX(envir)
  finalOutput = list(output = output, adj = prep$adj, orb = allGraphOrb)
  finalOutput
}

### This function prepares to optimise over the Ramsey polytope with specified bounds
prepareRamseyLP = function(numVertices, allGraphs, allGraphAdj, allGraphOrb, graphBounds, graphInfo, sparse = TRUE) {
  numVars = choose(numVertices, 2)
  auxMatrix = createPositionMatrix(numVertices = numVertices, symmetric = TRUE)
  graphBounds %<>%
    filter(!subsumed) %>%
    inner_join(graphInfo, by = join_by(graph == index, size == size)) %>%
    arrange(graph, direction) %>%
    mutate(numOpts = choose(numVertices, order), numCombinations = n_embed * numOpts, fullCombinations = numCombinations * size)
  numConstraints = sum(graphBounds$numCombinations)
  Dir = rep("", numConstraints)
  Rhs = rep(0,  numConstraints)
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
      if (is.null(allGraphOrb[[ind]])) {
        allGraphOrb[[ind]] = getOrbits(allGraphs[[ind]])
      }
      allGraphAdj[[ind]] = permuteGraph(allGraphs[[ind]], allGraphOrb[[ind]])
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
    curNEmbed = curInfo$n_embed[1]
    curNConst = curLength * curNOpts
    curNEntry = curNConst * curSize
    curNFull  = curLength * curCombos
    Rhs[pos + (1:curNFull)]   = rep(curInfo$bound,        curCombos)
    Dir[pos + (1:curNFull)]   = rep(curInfo$direction,    curCombos)
    miniIndices = rep(1:curNOpts, each = curLength)
    Map[pos + (1:curNFull), ] = cbind(rep(curInfo$number, curCombos), rep(1:curNEmbed, each = curNConst), rep(miniIndices, curNEmbed))
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
      allPermEdges = auxMatrix[cbind(as.vector(curPermEdges[1:curSize, ]), as.vector(curPermEdges[curSize + (1:curSize), ]))]
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
  output = list(mat = Mat, dir = Dir, rhs = Rhs, map = Map, adj = allGraphAdj)
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
### The return value is a 2m x p matrix (m: graph size, p: number of embeddings)
### The input is both the graph itself, and the output of getOrbits applied to it;
### however, if the latter is NULL, it is computed from scratch at the beginning.
permuteGraph = function(Graph, orbits = NULL) {
  Graph = matrix(Graph, ncol = 2)
  if (is.null(orbits)) {
    orbits = getOrbits(Graph)
  }
  n = orbits$n
  autoSize = orbits$autoSize
  if (autoSize > 1) {
    if (autoSize == factorial(n)) {
      allPerms = list(1:n)
    } else {
      goodComps = orbits$comps
      if (n <= MAX_SIZE) {
        allPerms = get(paste0("PERMS", n), envir = .GlobalEnv)
      } else {
        allPerms = permn(n)
      }
      for (ind in 1:nrow(goodComps)) {
        curComp = goodComps[ind, ]
        goodPerms = (map_int(allPerms, ~{.[curComp[1]]}) < map_int(allPerms, ~{.[curComp[2]]}))
        allPerms = allPerms[goodPerms]
      }
    }
  } else {
    allPerms = get(paste0("PERMS", n), envir = .GlobalEnv)
  }
  allPerms = matrix(unlist(allPerms), nrow = n)
  allPerms = rbind(allPerms[Graph[, 1], , drop = FALSE], allPerms[Graph[, 2], , drop = FALSE])
  allPerms
}

### This function computes a graph's full automorphism group from its generators,
### as well as all vertex and non-trivial possible edge orbits under its action.
### Note: for the edges, each orbit is listed in the format c(tailList, headList)
### The group is not returned, only the lexicographic comparisons it entails are.
getOrbits = function(Graph) {
  Graph  = matrix(Graph, ncol = 2)
  n = max(Graph)
  iGraph = graph_from_edgelist(Graph, directed = FALSE)
  autoGens = automorphism_group(iGraph)
  autoSize = as.integer(automorphisms(iGraph)$group_size)
  L = length(autoGens)
  if (L == 0) { return(list(orbitsV = split(1:n, 1:n), orbitsE = c(), comps = matrix(NA, 0, 2), autoSize = autoSize, n = n)) }
  if (autoSize == factorial(n)) { return(list(orbitsV = list(1:n), orbitsE = list(combn2(1:n)), comps = cbind(1:(n-1), 2:n), autoSize = autoSize, n = n)) }
  autoGens %<>% do.call(rbind, .)
  orbitsV = getMeet(autoGens)
  permGroup = constructPermGroup(autoGens, autoSize, returnID = FALSE)
  auxMatrix = createPositionMatrix(n, symmetric = TRUE)
  fullEdges = combn2(1:n)
  nC2 = nrow(fullEdges)
  allEdges = cbind(rep(1:nC2, autoSize - 1), auxMatrix[cbind(as.vector(permGroup[fullEdges[, 1], ]), as.vector(permGroup[fullEdges[, 2], ]))])
  auxG = graph_from_edgelist(allEdges, directed = FALSE)
  orbitsE = split(fullEdges, clusters(auxG)$membership)
  orbitsE = orbitsE[sapply(orbitsE, length) > 2]
  comps = getComparisons(permGroup)
  output = list(orbitsV = orbitsV, orbitsE = orbitsE, comps = comps, autoSize = autoSize, n = n)
  output
}

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

getComparisons = function(permGroup) {
  nr = nrow(permGroup)
  nc = ncol(permGroup)
  autoComp = apply(permGroup, 2, function(x) { min(which(x != 1:nr)) })
  autoNext = permGroup[cbind(autoComp, 1:nc)]
  goodComps = tibble(first = autoComp, second = autoNext) %>%
    distinct() %>%
    as.matrix()
  finalComps = getTransitiveReduction(goodComps)
  finalComps
}

### This function constructs a permutation group of given size from its generators
### The construction works by an implicit breadth-first search of the Cayley graph
### NOTE: The generators are specified by row, but the group is returned by column
constructPermGroup = function(generatorList, numElts, returnID = TRUE) {
  n = ncol(generatorList)
  L = nrow(generatorList)
  allElts = matrix(NA, numElts, n)
  allElts[1:(L + 1), ] = rbind(1:n, generatorList)
  curActive = 2:(L + 1)
  numActive = L
  pos = L + 1
  while (pos < numElts) {
    newActive = c()
    for (index in 1:L) {
      nextElts = allElts[curActive, generatorList[index, , drop = FALSE]]
      badInds = which(duplicated(rbind(allElts[1:pos, ], nextElts))) - pos
      if (length(badInds) < numActive) {
        newPos = pos + (1:(numActive - length(badInds)))
        allElts[newPos, ] = nextElts[setdiff(1:numActive, badInds), , drop = FALSE]
        pos = tail(newPos, 1)
        newActive = c(newActive, newPos)
      }
    }
    curActive = newActive
    numActive = length(curActive)
  }
  if (!returnID) { allElts = allElts[-1, , drop = FALSE] }
  allElts = t(allElts)
  allElts
}

### This function converts a permutation in array notation to a graph representation
permToGraph = function(perm, directed = TRUE) {
  extPerm = cbind(1:length(perm), perm)
  G = graph_from_edgelist(extPerm, directed = directed)
  G
}

### This function computes the meet of the cycle partitions of a set of permutations
getMeet = function(perms) {
  Graphs = apply(perms, 1, permToGraph)
  Union = do.call(igraph::union, Graphs)
  Parts = split(1:vcount(Union), clusters(Union)$membership)
  Parts
}

### This function creates a matrix M with M[i,j] being the edge position of (ij)
### The edges are ordered as follows: 12, 13, ..., 1N, 23, ..., 2N, ..., (N-1)N
### If symmetric = TRUE, the lower diagonal contains the same indices, reflected
createPositionMatrix = function(numVertices, symmetric = TRUE) {
  numPos = choose(numVertices, 2)
  allPos = combn2(1:numVertices)
  outMat = matrix(0, numVertices, numVertices)
  outMat[allPos] = 1:numPos
  if (symmetric) {
    outMat[allPos[, 2:1]] = 1:numPos
  }
  outMat
}

### This function tries to construct an optimized proof for an input bound (cut).
### LPfile contains the file from which it was proved; the graph is the one that
### is being bounded; sense is "min" for minimization or "max" for maximization.
### If eps is non-NA, the 1 (-1) increment in objective value is replaced by eps.
### The objective is to identify a proof with a minimal number of dual variables.
optimizeProof = function(LPfile, Graph, sense = "min", eps = NA) {
  envir = cplexAPI::openEnvCPLEX()
  model = cplexAPI::initProbCPLEX(envir, pname = "Ramsey Bound Optimisation")
  cplexAPI::readCopyProbCPLEX(envir, model, fname = LPfile)
  numVars  = cplexAPI::getNumCols(envir, model)
  numConst = cplexAPI::getNumRows(envir, model)
  numVerts = sqrt(2 * numVars + 1/4) + 1/2
  stopifnot(near(numVerts, round(numVerts)))
  numVerts = round(numVerts)
  n = max(Graph)
  stopifnot(n <= numVerts)
  auxMatrix = createPositionMatrix(numVerts, symmetric = TRUE)
  curObjective = rep(0, numVars)
  curObjective[curEdges] = 1
  cplexAPI::chgObjCPLEX(envir, model, ncols = numVars, ind = 0:(numVars - 1), val = curObjective)
  cplexAPI::setObjDirCPLEX(envir, model, ifelse(sense == "min", CPX_MIN, CPX_MAX))
  cplexAPI::setIntParmCPLEX(envir, cplexAPI::CPXPARAM_Preprocessing_Presolve, 0)
  cplexAPI::dualoptCPLEX(envir, model)
  curObjValue = cplexAPI::solutionCPLEX(envir, model)$objval
  objTarget = ifelse(sense == "min", ceiling(curObjValue), floor(curObjValue))
  dualLPFile = str_replace(LPfile, ".lp", "Dual.lp")
  dualWriteCPLEX(envir, model, fname = dualLPFile)
  cplexAPI::readProbCPLEX(envir, model, fname = dualLPfile)
  origDualObj = cplexAPI::getObjCPLEX(envir, model, 0, numConst - 1)
  epsilon = ifelse(is.na(eps), 1, eps)
  newBound = objTarget + ifelse(sense == "min", -1, 1) * (1 - epsilon)
  newSense = ifelse(sense == "min", "G", "L")
  newRange = rbind(1:numConst, numConst + (1:numConst), 1:numConst, numConst + (1:numConst))
  cplexAPI::newColsCPLEX(envir, model, ncols = numConst, lb = rep(0, numConst), xctype = rep(cplexAPI::CPX_CONTINUOUS, numConst))
  cplexAPI::addRowsCPLEX(envir, model, nrows = 2 * numConst + 1, nnz = 5 * numConst, rhs = c(rep(0, 2 * numConst), newBound), sense = c(rep("G", 2 * numConst), newSense), 
                        matbeg = c(rep(2, 2 * numConst), numConst), matval = c(as.vector(newRange), 1:numConst), matind = c(rep(c(-1, 1, 1, 1), numConst), origDualObj))
  newDualObj = rep(c(0, 1), each = numConst)
  cplexAPI::chgObjCPLEX(envir, model, ncols = numVars, ind = 0:(2 * numConst - 1), val = newDualObj)
  cplexAPI::setObjDirCPLEX(envir, model, ifelse(sense == "min", CPX_MAX, CPX_MIN))
  cplexAPI::primaloptCPLEX(envir, model)
  dualSolution = cplexAPI::solutionCPLEX(envir, model)$x[1:numConst]
  dualSolution
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
prepareBestProof = function(contradictions, listOfGraphs, labels = LETTERS, r = 3, s = 4, specialLabels = (if (r != s) { c("Z", "X") } else { "X" }), 
                            fname = paste0("ShortProofR", r, "S", s, ".RData"), plotGraphs = FALSE) {
  endpoint = 1 + (r != s)
  numGraphs = sapply(contradictions, function(x) { n_distinct(x$graph) })
  numBounds = sapply(contradictions, nrow)
  numSteps  = sapply(contradictions, function(x) { n_distinct(x$iter) })
  contraTab = tibble(index = 1:length(contradictions), use = numGraphs, bounds = numBounds, steps = numSteps) %>%
    arrange(use, bounds, steps)
  bestProof = contraTab$index[1]
  bestContradiction = contradictions[[bestProof]]
  bestSize = numGraphs[bestProof]
  stopifnot(length(specialLabels) == endpoint)
  usableLabels = setdiff(labels, specialLabels)
  if (bestSize > length(usableLabels) + endpoint) { print("Not enough labels!"); return() }
  uGraphs = sort(unique(bestContradiction$graph))
  firstGraphs = bestContradiction %>%
    slice(1:2) %>%
    pull(graph)
  nameTab = tibble(index = setdiff(uGraphs, firstGraphs), label = usableLabels[1:(bestSize - endpoint)]) %>%
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
