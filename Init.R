source("Ramsey.R")

# U = iterateLPProof(n = 9L, r = 3L, s = 4L, extremeOnly = TRUE)
# save(U, file = "ProofR4S3ExtremeOnly.RData")
# Qu = findMinimalGraphSet(extremeOnly = TRUE, inds = c(4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 27, 29, 30, 33, 35, 36, 37))
# save(Qu, file = "ResultsR3S4MinimalExtremeOnly.RData")
# load("ResultsR3S4MinimalExtremeOnly.RData")
# ResU = optimiseProofs(Qu[[2]])
# save(ResU, file = "OptimisedProofsR3S4MinimalExtreme.RData")

# Qv = findMinimalGraphSet(extremeOnly = TRUE, inds = c(5, 10, 11, 15, 22, 24, 29, 30)) # The original set leading to contradiction being lost we simply add 30 at the end
# save(Qv, file = "ResultsR3S4MinimalExtremeOnlyInitial.RData")
# load("ResultsR3S4MinimalExtremeOnlyInitial.RData")
# ResV = optimiseProofs(Qv[[2]])
# save(ResV, file = "OptimisedProofsR3S4MinimalExtremeInitial.RData")

# W = iterateLPProof(n = 9L, r = 3L, s = 4L, maxOrder = 8L, minAuto = 24L, lowerOnly = TRUE, strongCuts = TRUE)
# save(W, file = "ProofR4S3Auto24LowerOnly.RData")
# Qw = findMinimalGraphSet(maxOrder = 8L, minAuto = 24L, lowerOnly = TRUE, inds = c(3,4,5,6,7,8,9,10,12,13,14,15,25,28,34,35,37,38,39,42,44,48,49,50,52,55,56,58,60,
#                                                                                    61,65,67,69,70,74,75,79,81,84,89,90,94,96,98,100,103,107,110,114,115,116,118,119,
#                                                                                    120,122,126,128,129,130,133,134,137,144,149,150,151,152,154,155,156,162,163,165,168,
#                                                                                    172,174,176,183,184,190,194,198,202,204,217,218,220,221,223,226,227,236,237,240,242,
#                                                                                    246,248,250,253,258,266,283,289,290,300,302,310,312,318,323,334,338,341))
# save(Qw, file = "ResultsR3S4MinimalAuto24LowerOnly.RData")
# ResW = optimiseProofs(Qw[[2]])
# save(ResW, file = "OptimisedProofsR3S4MinimalAuto24LowerOnly.RData")

# X = iterateLPProof(n = 9L, r = 3L, s = 4L, maxOrder = 8L, treesOnly = TRUE)
# save(X, file = "ProofR4S3Trees.RData")
# Qx = findMinimalGraphSet(treesOnly = TRUE,  maxOrder = 8L, inds = c(4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20, 22, 23, 24, 27, 32, 34, 35, 38, 42, 43, 46, 49))
# save(Qx, file = "ResultsR3S4MinimalTreesUpTo8.RData")
# load("ResultsR3S4MinimalTreesUpTo8.RData")
# ResX = optimiseProofs(Qx[[2]])
# save(ResX, file = "OptimisedProofsR3S4MinimalTreesUpTo8.RData")

# Y = iterateLPProof(n = 9L, r = 3L, s = 4L, maxOrder = 6L, lowerOnly = TRUE)
# save(Y, file = "ProofR4S3LowerOnly.RData")
# Qy = findMinimalGraphSet(maxOrder = 6L, lowerOnly = TRUE, inds = c(4,5,7,8,9,10,11,13,14,15,18,19,20,21,23,24,25,27,28,29,30,31,32,33,35,36,41,42,47,48,50,51,
#                                                                   54,55,61,63,66,67,72,73,77,81,87,88,90,91,95,96,99,100,104,106,114,116,117,118,119,120,122,
#                                                                   123,125,129,130,133,134,138,139,141,142))
# save(Qy, file = "ResultsR3S4MinimalLowerOnly.RData")
# load("ResultsR3S4MinimalLowerOnly.RData")
# ResY = optimiseProofs(Qy[[2]])
# save(ResY, file = "OptimisedProofsR3S4MinimalLowerOnly.RData")

# Z = iterateLPProof(n = 9L, r = 3L, s = 4L, maxOrder = 6L, partiteOnly = TRUE)
# save(Z, file = "ProofR4S3PartiteOnly.RData")
# Qz = findMinimalGraphSet(maxOrder = 6L, partiteOnly = TRUE, inds = c(4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 27))
# save(Qz, file = "ResultsR3S4MinimalBipartiteOnly.RData")
# load("ResultsR3S4MinimalBipartiteOnly.RData")
# ResZ = optimiseProofs(Qz[[2]])
# save(ResZ, file = "OptimisedProofsR3S4MinimalBipartiteOnly.RData")

# Z = iterateLPProof(n = 13L, r = 3L, s = 5L, maxOrder = 12L, treesOnly = TRUE)
# save(Z, file = "ProofR5S3TreesOnly.RData")
# Z = iterateLPProof(n = 18L, r = 4L, s = 4L, maxOrder = 12L, twoTreesOnly = TRUE)
# save(Z, file = "ProofR4S4TwoTreesOnly.RData")

# noMult3 = iterateLPProof(n = 6L, r = 3L, s = 3L, maxOrder = 5L, inds = c(2, 3, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 30))
# mult6 = iterateLPProof(n = 6L, r = 3L, s = 3L, maxOrder = 5L, minAuto = 6, factor = TRUE)

# T = iterateLPProof(n = 9L, r = 3L, s = 4L, maxOrder = 8L, minAuto = 24L, factor = TRUE)
# save(T, file = "ProofR4S3MinimalAuto24Factor.RData")
Qt = findMinimalGraphSet(maxOrder = 8L, minAuto = 24L, factor = TRUE, inds = c(3,4,5,7,11,12,13,20,23,35,37,39,41,42,44,46,48,49,50,51,52,59,68,69,82,86,93,95,101,102,
                                                                               104,107,109,114,119,121,122,123,124,128,139,143,145,152,159,166,170,171,181,189,197,198,
                                                                               204,209,210,216,219,223,227,231,236,240,247,253,257,261,262,266,270,275,276,282,284,286,
                                                                               288,291,292,294,296,299,302))
save(Qt, file = "ResultsR3S4MinimalAuto24Factor.RData")