#' DTM Triplet Statistic for Cross-group Testing of Mean Proportions
#'
#' \code{nodetest} returns Dirichlet-tree multinomial (DTM) node and triplet test 
#' statistics based on Method-of-Moments (MoM) asymptotic test on each internal node.
#' 
#' @param pstrct An object returned from function \code{phylostructure}.
#' @param group.data A list where each element is a matrix of taxonomic 
#' counts(columns) for each sample(rows). The ordering of columns needs to be the 
#' same as the ordering of leaves in the phylogenetic tree.
#' @return A list containing the following components:
#' \item{\code{MoMp}}{MoM p-value for each internal node
#' using the DTM model.}
#' \item{\code{tripletstat}}{Triplet test statistics under the DTM model. Each 
#' row corresponds to a triplet, where the first three columns contain the node 
#' labels and the last column contains the triplet statistic.}
#' \item{\code{w}}{Maximum of all triplet statistics.}
#' @author Yunfan Tang
#' @references Tang, Y., Ma, Li. and Nicolae, D. L. (2017). Phylogenetic 
#' Dirichlet-multinomial model for microbiome data. 
#' \href{https://arxiv.org/abs/1610.08974}{arXiv:1610.08974} [stat.AP].
#' @references Wang, T. and Zhao, H. (2017). A Dirichlet-tree multinomial regression 
#' model for associating dietary nutrients with gut microorganisms. Biometrics.
#' @examples 
#' library(ape)
#' set.seed(10)
#' 
#' ## Generate a random binary phylogenetic tree with 10 leaf nodes
#' pstrct <- phylostructure(rtree(10))
#' 
#' ## Simulate microbiome samples for two groups from multinomial distribution
#' p1 <- c(rep(0.12, 3), rep(0.08, 3), rep(0.1, 4))
#' p2 <- p1 + 0.001 * c(c(1, -1), rep(0, 8))
#' n <- 1000 #Number of sequences in each sample
#' m <- 200 #Number of samples in each group
#' group.data <- list(x1 = t(rmultinom(m, n, p1)), x2 = t(rmultinom(m, n, p2)))
#' 
#' ## Calculte triplet statistics
#' nodetest(pstrct, group.data)
#' 
#' @importFrom HMP Xmcupo.sevsample
#' @importFrom stats dchisq integrate pchisq qchisq
#' @export

nodetest <- function(pstrct, group.data) {
  phylotree <- pstrct$phylotree
  phyloparent <- pstrct$phyloparent
  phylochildren <- pstrct$phylochildren
  
  MoMp <- matrix(0, phylotree$Nnode, 2, dimnames = list(c(), c("Node", "p-value")))
  for (gd in group.data) if (ncol(gd) != ntaxa(pstrct$phylotree)) 
    stop("Number of columns in group.data is not equal to number of leaves")
  for (i in 1:phylotree$Nnode) {
    NodeID <- i + ntaxa(phylotree)
    c1 <- phylochildren[NodeID, 1]
    c2 <- phylochildren[NodeID, 2]
    gp <- lapply(group.data, function(g) {
      k1 <- g[, pstrct$descendant[c1, ]]
      k2 <- g[, pstrct$descendant[c2, ]]
      if (is.matrix(k1)) 
        k1 <- rowSums(k1)
      if (is.matrix(k2)) 
        k2 <- rowSums(k2)
      tm <- cbind(k1, k2)
      tm[rowSums(tm) != 0, ]
    })
    MoMp[i, ] <- c(NodeID, Xmcupo.sevsample(gp)$`p value`)
  }
  
  ts <- matrix(0, phylotree$Nnode, 4, dimnames = list(
    c(), c("Node1", "Node2", "Node3", "stat")))
  for (i in 2:phylotree$Nnode) {
    NodeID <- i + ntaxa(phylotree)
    Par <- phyloparent[NodeID]
    GPar <- phyloparent[Par]
    if (GPar == 0) 
      next
    ts[i, 1:3] <- c(GPar, Par, NodeID)
    ts[i, 4] <- sum(qchisq(1 - MoMp[c(NodeID, Par, GPar) - ntaxa(phylotree), 2], 1))
  }
  ts <- ts[rowSums(ts[, 1:3]) != 0, ]
  list(MoMp = MoMp, tripletstat = ts, w = max(ts[, 4]))
}
