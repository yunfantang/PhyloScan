#' Generating Node Relations
#'
#' \code{phylostructure} returns the neighboring nodes and descendant leaves 
#' for each internal node of a binary phylogenetic tree.
#' @param phylotree A \code{phylo-class} object as in the \code{phyloseq} package. 
#' The tree has to be binary.
#' @details For each internal node (i.e. non-leaf node) in \code{phylotree}, 
#' this function calculates its parent node, child nodes and all of its descendant
#' leaves. The labeling of nodes are the same as in \code{phylotree$edge}, which is 
#' assumed to obey the following convention. Leaves are labeled from \code{1} to 
#' \code{ntaxa(phylotree)} and internal nodes are labeled from 
#' \code{ntaxa(phylotree) + 1} to \code{ntaxa(phylotree) + phylotree$Nnode}. The 
#' root label is \code{ntaxa(phylotree) + 1}. Parent node always has a lower label 
#' than its child internal nodes.
#' @return A list containing the following components:
#' \item{\code{phylotree}}{Input object}
#' \item{\code{phylochildren}}{Numerical matrix of child node labels. The kth row 
#' contains the children node labels of node k. Leaf nodes have their their row 
#' vectors set to zero.}
#' \item{\code{phyloparent}}{Numerical vector of parent node labels. The kth element 
#' contains the parent node label of node k. Root node has its parent set to zero.}
#' \item{\code{descendant}}{Logical matrix of decendant leaves. The kth row 
#' shows whether the leaves are descendants of node k.}
#' @author Yunfan Tang
#' @examples
#' library(ape)
#' set.seed(10)
#' 
#' ## Analyze a random binary phylogenetic tree with 6 leaf nodes
#' pstrct <- phylostructure(rtree(6))
#' ## Show children of node 7 (root):
#' pstrct$phylochildren[7, ]
#' ## Show parent of node 1:
#' pstrct$phyloparent[1]
#' ## List all leaf nodes under node 8:
#' which(pstrct$descendant[8, ])
#' 
#' @export

phylostructure <- function(phylotree) {
  if(class(phylotree) != "phylo")
    stop("phylotree is not a phylo class")
  phyloparent <- numeric(phylotree$Nnode + ntaxa(phylotree))
  phylochildren <- matrix(0, phylotree$Nnode + ntaxa(phylotree), 2)
  for (i in 1:nrow(phylotree$edge)) {
    i1 <- phylotree$edge[i, 1]
    i2 <- phylotree$edge[i, 2]
    if (i1 <= ntaxa(phylotree)) 
      stop("Internal node label is not larger than ntaxa(phylotree)")
    if (i2 > ntaxa(phylotree) && i1 > i2) 
      stop("Parent node label is larger than child internal node label")
    if (all(phylochildren[i1, ] > 0)) 
      stop("Phylogenetic tree is not binary")
    
    phyloparent[i2] <- i1
    if (phylochildren[i1, 1] == 0) 
      phylochildren[i1, 1] <- i2 else phylochildren[i1, 2] <- i2
  }
  
  descendant <- matrix(FALSE, phylotree$Nnode + ntaxa(phylotree), ntaxa(phylotree))
  for (i in 1:ntaxa(phylotree)) descendant[i, i] <- TRUE
  processed <- logical(phylotree$Nnode + ntaxa(phylotree))
  processed[1:ntaxa(phylotree)] <- TRUE
  while (!all(processed)) {
    for (i in (ntaxa(phylotree) + 1):(ntaxa(phylotree) + phylotree$Nnode)) {
      if (all(processed[phylochildren[i, ]])) {
        descendant[i, descendant[phylochildren[i, 1], ]] <- TRUE
        descendant[i, descendant[phylochildren[i, 2], ]] <- TRUE
        processed[i] <- TRUE
      }
    }
  }
  list(phylotree = phylotree, phylochildren = phylochildren, phyloparent = phyloparent, 
       descendant = descendant)
}
