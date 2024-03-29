MAX_SIZE  = 9L
print(paste("Initializing permutations up to size", MAX_SIZE))
for (ind in 1:MAX_SIZE) {
  print(ind)
  assign(paste0("PERMS", ind), permn(ind), envir = .GlobalEnv)
}

findAllUE = function(minOrder = 3L, maxOrder = 8L) {
  goodPairs = matrix(NA, 10^(maxOrder - minOrder + 1L), 2)
  pos = 1
  allGraphs = list()
  for (gorder in minOrder:maxOrder) { 
    allGraphs %<>% 
      c(getGraphs(gorder)) 
  }
  graphTab = tibble(index = 1:length(allGraphs), size = map_int(allGraphs, length)/2, redSize = size - 2L, order = map_int(allGraphs, max), 
                    minDeg = map_int(allGraphs, ~{min(table(.))}), maxDeg = map_int(allGraphs, ~{max(table(.))}))
  allComp = inner_join(graphTab, graphTab, by = join_by(size <= redSize, order <= order, minDeg <= minDeg, maxDeg <= maxDeg)) %>%
    select(index.x, index.y, size.x, size.y)
  smallPrimes = c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L, 41L, 43L, 47L, 53L, 59L, 61L, 67L, 71L, 73L)
  maxSize = max(allComp$size.y)
  for (prime in smallPrimes) {
    if (prime <= maxSize) {
      allComp %<>% 
        mutate_at("rp", ~{. & (mod(size.x, prime) != 0 | mod(size.y, prime) != 0)})
      print(prime)
    }
  }
  allComp %<>%
    filter(rp)
  allGraphs %<>%
    sapply(makeIgraph)
  numComp = nrow(allComp)
  for (index in 1:numComp) {
    if (index %% 10000 == 0) { print(index) }
    curPair = allComp %>%
      slice(index)
    ind1 = curPair$index.x
    ind2 = curPair$index.y
    if (checkUE(allGraphs[[ind1]], allGraphs[[ind2]])) {
      goodPairs[pos,] = c(ind1, ind2)
      pos = pos + 1
    }
  }
  goodPairs
}

### This auxiliary function checks if the subgraph is uniformly embeddable into
### the graph, meaning that every edge is represented the same number of times
checkUE = function(Graph, Subgraph) {
  if ((vcount(Subgraph) > vcount(Graph)) || (vcount(Subgraph) == vcount(Graph) && ecount(Subgraph) >= ecount(Graph))) {
    return(FALSE)
  }
  check = subgraph_isomorphisms(pattern = Subgraph, target = Graph)
  L = length(check)
  m = ecount(Subgraph)
  M = ecount(Graph)
  N = vcount(Graph)
  estK = (L * m) / M
  if (L == 0 || !near(estK, round(estK))) { 
    return(FALSE)
  }
  subEdges = as_edgelist(Subgraph)
  fullMap = lapply(check, function(x) { matrix(x[subEdges], ncol = 2) }) %>%
    do.call(rbind, .)
  swapOrder = which(fullMap[,1] > fullMap[,2])
  if (length(swapOrder) > 0) {
    fullMap[swapOrder, ] = fullMap[swapOrder, 2:1]
  }
  sparseMat = sparseMatrix(i = fullMap[, 1], j = fullMap[, 2], x = 1, dims = c(N, N), index1 = TRUE, use.last.ij = FALSE)
  result = all(near(sparseMat@x, estK))
  result
}

### This auxiliary function converts an igraph into a flat graph representation
flattenGraph = function(iGraph) {
  eList = get.edgelist(iGraph)
  output = as.vector(eList[order(eList[, 1], eList[, 2]), ])
  output
}

### This auxiliary function converts a flat graph representation into an igraph
makeIgraph = function(Graph) {
  output = graph_from_edgelist(matrix(Graph, ncol = 2), directed = FALSE)
  output
}

### This auxiliary function converts an igraph into its canonical representative
canonicallyOrderGraph = function(iGraph) {
  output = permute(iGraph, canonical_permutation(iGraph)$labeling)
  output
}

### This auxiliary function makes all graphs in a given input file canonical
canonicallyOrderAllGraphs = function(graphFile, varName = "allG") {
  load(graphFile)
  allGraphs = get(varName)
  L = length(allGraphs)
  print(paste("There are", L, "graphs to process"))
  for (ind in 1:L) {
    if (ind %% 10000 == 0) { print(ind) }
    curGraph = allGraphs[[ind]] %>%
      makeIgraph %>%
    canonicallyOrderGraph %>%
    flattenGraph
    allGraphs[[ind]] = curGraph
  }
  assign(varName, allGraphs)
  save(list = varName, file = graphFile)
  allGraphs
}

### This function extracts all the constrained graphs from a given list of files
### varName specifies the name of the variable under which the results are saved
### These results are expected to be the outputs of the iterateLPProof function
combineUsefulGraphs = function(listOfFiles = paste0("ProofR5S3", c("ET", "Trees"), "OnlyPart2NoContradiction.RData"), varName = "Z", outFile = "Custom.RData") {
  Len = length(listOfFiles)
  allRes = vector("list", Len)
  for (ind in 1:Len) {
    print(ind)
    fname = listOfFiles[[ind]]
    load(fname)
    allRes[[ind]] = get(varName)
    rm(list = varName)
  }
  usefulGraphs = vector("list", Len)
  for (ind in 1:Len) {
    curRes = allRes[[ind]]
    allGraphs = curRes$graphs
    allBounds = curRes$boundTab
    graphInds = sort(unique(allBounds$graph))
    useGraphs = allGraphs[graphInds]
    usefulGraphs[[ind]] = useGraphs
  }
  usefulGraphs %<>% 
    do.call(c, .)
  badInds = which(duplicated(usefulGraphs))
  allG = usefulGraphs[-badInds]
  save(allG, file = outFile)
  allG
}

### This function starts with a density at startN and computes the corresponding
### density up to stopN; rounding occurs at each step to account for integrality
### If lower = TRUE, this is a lower density, otherwise, it is an upper density.
### The density should be specified by 2 integers, ie. (numerator, denominator).
densityRatchet = function(density, startN, stopN, lower = TRUE) {
  stopifnot(stopN >= startN)
  stopifnot(all(density >= 0) && density[1] <= density[2])
  dMatrix = tibble(ind = startN:stopN, num = density[1], den = density[2]) %>%
    mutate_all(as.integer)
  pos = startN
  while (pos <= stopN) {
    density[1] %<>% multiply_by(choose(pos, 2))
    newNum = divide_by_int(density[1], density[2])
    if (lower && (density[1] %% density[2] != 0)) { 
      newNum %<>% add(1L)
    }
    density = c(newNum, choose(pos, 2))
    GCD = pracma::gcd(density[1], density[2])
    density %<>% divide_by(GCD)
    dMatrix[pos - startN + 1, ] = tibble(ind = pos, num = density[1], den = density[2])
    pos = pos + 1
  }
  dMatrix %<>% 
    mutate(ratio = num/den)
  dMatrix
}

### This function starts with a target density at stopN and finds the required
### density at startN; rounding occurs at each step to account for integrality
### If lower = TRUE, this is a lower density, otherwise, it is an upper density.
### The target should be specified by 2 integers, ie. (numerator, denominator).
inverseRatchet = function(target, startN, stopN, lower = TRUE) {
  stopifnot(stopN >= startN)
  stopifnot(all(target >= 0) && target[1] <= target[2])
  dMatrix = tibble(ind = startN:stopN, num = target[1], den = target[2]) %>%
    mutate_all(as.integer)
  pos = stopN
  while (pos > startN) {
    target[1] %<>% multiply_by(choose(pos, 2))
    newNum = divide_by_int(target[1], target[2])
    if (target[1] %% target[2] == 0) { 
      newNum %<>% add(ifelse(lower, -1L, 1L))
    }
    target = c(newNum, choose(pos, 2))
    GCD = pracma::gcd(target[1], target[2])
    target %<>% divide_by(GCD)
    dMatrix[pos - startN, ] = tibble(ind = pos - 1, num = target[1], den = target[2])
    pos = pos - 1
  }
  dMatrix %<>% 
    mutate(ratio = num/den)
  dMatrix
}

### This function extracts connected graphs from the files made by McKay's nauty
### Make sure that the conversion is done with the -el0o1 option in nauty::showg
### If extremeOnly = TRUE, only keeps the graphs with sizes within 1 of min/max.
getGraphs = function(numVerts = 2, extremeOnly = FALSE, treesOnly = FALSE, twoTreesOnly = FALSE, ETOnly = FALSE, partiteOnly = FALSE) {
  if (!treesOnly && !twoTreesOnly && !ETOnly && !partiteOnly) {
    fn = paste0("Graphs", numVerts, ".RData")
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
    curText = readLines(paste0(ifelse(treesOnly, "trees", ifelse(twoTreesOnly, "TwoTrees", ifelse(ETOnly, "ET", ifelse(partiteOnly, "Bipartite", "graph")))), numVerts, "c", ".txt"))
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
    if (extremeOnly) {
      curCounts = sapply(allG, length)/2
      goodPos = ((curCounts <= numVerts) | (curCounts >= choose(numVerts, 2) - 1))
      allG = allG[goodPos]
    }
  }
  allG
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
    curGraph    = makeIgraph(allGraphs[[curRow$graph]])
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
    subsumedBySuper = sapply(allGraphs[subsumedSuperCandidates$graph], function(x) { subgraph_isomorphic(curGraph, makeIgraph(x)) })
    subsumedBySub   = sapply(allGraphs[  subsumedSubCandidates$graph], function(x) { subgraph_isomorphic(makeIgraph(x), curGraph) })
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
  goodRows   = which(!near(dualVariables, 0, tol = EPSILON))
  dualValues = dualVariables[goodRows]
  dualValues = dualValues/min(abs(dualValues))
  dir        = LP$dir[goodRows]
  relInfo    = LP$map[goodRows, , drop = FALSE]
  if (shortProofs) {
    proof     = tibble(coeff = dualValues, constr = goodRows, dir = dir, bound = relInfo[,1], version = relInfo[,2], combo = relInfo[,3])
  } else {
    subMatrix = as.matrix(LP$mat[goodRows, , drop = FALSE])
    rhs       = LP$rhs[goodRows]
    proof     = cbind(coeff = dualValues, subMatrix, rhs = rhs)
    flipRows  = which(dir == ifelse(lower, "L", "G"))
    proof[flipRows, ] %<>% magrittr::multiply_by(-1)
    proof %<>% cbind(bound = relInfo[,1])
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
    if (n <= MAX_SIZE) {
      allPerms = get(paste0("PERMS", n), envir = .GlobalEnv)
    } else {
      allPerms = permn(n)
    }
  }
  allPerms = matrix(unlist(allPerms), nrow = n)
  allPerms = rbind(allPerms[Graph[, 1], , drop = FALSE], allPerms[Graph[, 2], , drop = FALSE])
  allPerms
}

### This function computes a graph's full automorphism group from its generators,
### as well as all vertex and non-trivial possible edge orbits under its action.
### Note: the edge orbits are numbered with the universal edge-numbering scheme.
### The group is not returned, only the lexicographic comparisons it entails are.
getOrbits = function(Graph) {
  n = max(Graph)
  iGraph = makeIgraph(Graph)
  nC2 = choose(n, 2)
  autoGens = automorphism_group(iGraph)
  autoSize = as.integer(automorphisms(iGraph)$group_size)
  L = length(autoGens)
  if (L == 0) { return(list(orbitsV = split(1:n, 1:n), orbitsE = c(), comps = matrix(NA, 0, 2), autoSize = autoSize, n = n)) }
  fullEdges = createPositionMatrix(n, twoColumn = TRUE)
  if (autoSize == factorial(n)) { return(list(orbitsV = list(1:n), orbitsE = list(fullEdges), comps = cbind(1:(n-1), 2:n), autoSize = autoSize, n = n)) }
  autoGens %<>% do.call(rbind, .)
  orbitsV = getMeet(autoGens)
  auxMatrix = createPositionMatrix(n, symmetric = TRUE, twoColumn = FALSE)
  allEdges = cbind(rep(1:nC2, each = L), auxMatrix[cbind(as.vector(autoGens[, fullEdges[, 1], drop = FALSE]), as.vector(autoGens[, fullEdges[, 2], drop = FALSE]))])
  auxG = makeIgraph(allEdges)
  orbitsE = split(fullEdges, clusters(auxG)$membership)
  orbitsE = orbitsE[sapply(orbitsE, length) > 2]
  permGroup = constructPermGroup(autoGens, autoSize)
  comps = getComparisons(permGroup)
  output = list(orbitsV = orbitsV, orbitsE = orbitsE, comps = comps, autoSize = autoSize, n = n)
  output
}

### This function constructs a minimal set of constraints of the form s[a] < s[b]
### that the lexicographically smallest element of a coset of its input satisfies
### The input is the permutation group with one column per element minus identity
### The output is a two-column matrix where a row [ab] corresponds to s[a] < s[b]
### NOTE: assumes that if the identity is included, then it has to be in column 1
getComparisons = function(permGroup) {
  nr = nrow(permGroup)
  nc = ncol(permGroup)
  if (all(permGroup[, 1] == 1:nr)) {
    permGroup = permGroup[, -1, drop = FALSE]
    nc = nc - 1
  }
  autoComp = apply(permGroup, 2, function(x) { min(which(x != 1:nr)) })
  autoNext = permGroup[cbind(autoComp, 1:nc)]
  goodComps = tibble(first = autoComp, second = autoNext) %>%
    distinct()
  finalComps = goodComps %>%
    group_by(second) %>%
    filter(first == max(first)) %>%
    ungroup() %>%
    arrange(first, second) %>%
    as.matrix()
  finalComps
}

### This function constructs a permutation group of given size from its generators
### The construction works by an implicit breadth-first search of the Cayley graph
### NOTE: The generators are specified by row, but the group is returned by column
constructPermGroup = function(generatorList, numElts = NA) {
  n = ncol(generatorList)
  L = nrow(generatorList)
  if (is.na(numElts)) { numElts = factorial(n) }
  allElts = matrix(NA, numElts, n)
  allElts[1:(L + 1), ] = rbind(1:n, generatorList)
  curActive = 2:(L + 1)
  numActive = L
  pos = L + 1
  while (numActive > 0) {
    newActive = c()
    for (index in 1:L) {
      nextElts = allElts[curActive, generatorList[index, , drop = FALSE], drop = FALSE]
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
  allElts = allElts[1:pos, , drop = FALSE]
  allElts = t(allElts)
  allElts
}

### This function converts a permutation in array notation to a graph representation
permToGraph = function(perm) {
  extPerm = cbind(1:length(perm), perm)
  G = makeIgraph(extPerm)
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
### The edges are ordered as: 12, 13, 23, 14, 24, 34, ..., 1N, 2N, ..., (N-1)N
### If symmetric = TRUE, the lower diagonal contains the same indices, reflected
### If twoColumn = TRUE, returns the list of edges in two-column format instead!
createPositionMatrix = function(numVertices, symmetric = TRUE, twoColumn = FALSE) {
  numPos = choose(numVertices, 2)
  allPos = (numVertices + 1) - combn2(1:numVertices)
  allPos = allPos[numPos:1, 2:1]
  if (twoColumn) {
    return(allPos)
  }
  outMat = matrix(0, numVertices, numVertices)
  outMat[allPos] = 1:numPos
  if (symmetric) {
    outMat[allPos[, 2:1]] = 1:numPos
  }
  outMat
}

### This function returns the value of N > 0 such that numPairs equals N choose 2
### Returns an error if the resulting solution is not an integer within tolerance
invertNumPairs = function(numPairs) {
  baseNum = (sqrt(numPairs * 8 + 1) + 1)/2
  baseInt = as.integer(baseNum)
  stopifnot(near(baseNum, baseInt, tol = EPSILON))
  baseInt
}

### This function traces back through a list of short proofs to obtain contradictions
reconstructContradictions = function(allProofs, boundTab, n, specialIndex = NULL) {
  lastIter = max(boundTab$iter)
  possibleEdges = choose(n, 2)
  lastProofs = tibble()
  if (!is.null(specialIndex)) {
    specialBound = boundTab %>% 
      filter(graph == specialIndex & direction == "G") %>%
      mutate(density = bound/size) %>%
      mutate(expectedSize = density * possibleEdges) %>%
      select(number, direction, expectedSize)
    lastProofs = inner_join(specialBound, specialBound, by = join_by(expectedSize >= expectedSize))
  }
  extraProofs = boundTab %>%
    filter(iter %in% c(1, lastIter)) %>%
    mutate(density = bound/size) %>%
    mutate(expectedSize = density * possibleEdges) %>%
    select(number, direction, expectedSize)
  extraProofsG = extraProofs %>%
    filter(direction == "G") %>%
    mutate_at("expectedSize", ceiling)
  extraProofsL = extraProofs %>%
    filter(direction == "L") %>%
    mutate_at("expectedSize", floor)
  lastProofs = bind_rows(lastProofs, inner_join(extraProofsG, extraProofsL, by = join_by(expectedSize > expectedSize)))
  numContradictions = nrow(lastProofs)
  allContradictions = vector("list", numContradictions)
  for (index in 1:numContradictions) {
    curRow = lastProofs %>%
      slice(index)
    curIndices = unique(c(curRow$number.x, curRow$number.y))
    allIndices = curIndices
    contradiction = vector("list", lastIter)
    pos = 1
    ### Carrying out a breadth-first-like search to reach a contradiction
    while (any(curIndices > 2)) {
      contradiction[[pos]] = boundTab %>% 
        slice(curIndices)
      pos = pos + 1
      curIndices = setdiff(curIndices, 1:2)
      if (length(curIndices) > 0) {
        curProofs  = allProofs[curIndices]
        curIndices = lapply(curProofs, function(x) {x$bound}) %>%
          unlist() %>%
          c() %>%
          unique()
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
### If plotGraphs = TRUE, the graphs are also plotted into PDF files, each named as [label].pdf
prepareBestProof = function(contradictions, listOfGraphs, labels = outer(LETTERS, letters, paste0), r = 3, s = 4, specialLabels = (if (r != s) { c("Zz", "Xx") } else { "Xx" }), 
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
    slice(1:endpoint) %>%
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
    allGraphObj = sapply(listOfGraphs, makeIgraph)
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
  if (!is.null(fname)) {
    save(bestContradiction, miniContradiction, file = fname)
  }
  bestContradiction
}

### This auxiliary function fills a matrix with zero rows up to a specified size
fillWithZeros = function(x, size) {
  Mat = matrix(0, size, ncol(x))
  Mat[1:nrow(x), ] = x
  Mat
}
