#' Map Mass2Motifs on mass spectral molecular network
#'
#' @param edges edges file from GNPS
#' @param motifs motif summary table from MS2LDA
#' @param prob minimal probability score for a Mass2Motif to be included. Default is 0.01. 
#' @param overlap minimal overlap score for a Mass2Motif to be included. Default is 0.3.
#' @param top parameter specifiying how many most shared motifs per molecular family (network component index) should be shown. Default is 5.
#'
#' @return edges and nodes file with Mass2Motifs mapped
#' @export
#'
#' @examples
#' @import stats 

Mass2Motif_2_Network <- function(edges,motifs,prob = 0.01,overlap = 0.3, top = 5){
  
  if (colnames(motifs)[1] != "scans"){
    print("WARNING: First column is used for ID matching")
    colnames(motifs)[1] <- "scans"
  }
   
  # set cutoff for motifs to be included: Probability min. 0.01 and Overlap min 0.3 is default
  motifs <- motifs[intersect(which(motifs$probability >= prob), which(motifs$overlap >= overlap)),]
  
  # create additional column in edges file containing shared motifs between each node pair
  shared_motifs <- function(nodes){
    a <- motifs$motif[motifs$scans %in% nodes[1]]
    b <- motifs$motif[motifs$scans %in% nodes[2]]
    out <- paste(sort(intersect(a,b)),collapse = ",")
    return(out)
  }
  
  edges$SharedMotifs <- apply(edges[,c(1:2)], 1, shared_motifs)
  edges$interact <- "cosine"
  
  # add additional rows for each shared motif, so each shared motif can be displayed with an individual edge
  l <- strsplit(edges$SharedMotifs,split=",")
  
  edges_m <- edges[rep(seq_len(dim(edges)[1]), lengths(l)), ]
  edges_m$interact <- unlist(l)
  edges <- rbind(edges, edges_m)
  
  # add additional column in edges file containing the x most shared motifs per molecular family
  agg <- stats::aggregate(interact~ComponentIndex, data = edges_m[-which(edges_m$ComponentIndex == -1),], paste0, collapse=",")
  
  agg_c <- strsplit(as.character(agg$interact),split = ",")
  c <- lapply(agg_c,plyr::count)
  topX <- lapply(c, function(x) x[order(x$freq,decreasing=T), ])
  topX <- lapply(topX, function(x) x[1:top,1])
  topX <- unlist(lapply(topX, paste0, collapse = ","))
  agg$topX <- topX
  
  edges$topX <- agg$topX[match(edges$ComponentIndex,agg$ComponentIndex)]
  
  # reorder columns
  edges <- edges[,c("CLUSTERID1", "interact", "CLUSTERID2", "DeltaMZ", "MEH", "Cosine", 
                                  "OtherScore", "ComponentIndex", "SharedMotifs", "topX")]
  
  edges <- edges[order(edges$ComponentIndex),]
  
  # create node table containing overlap scores of motifs per node
  motifs[-1] = apply(motifs[-1],2,as.character)
  motifs_cytoscape <- stats::aggregate(motifs[-1],by=list(motifs$scans),c)
  
  ul <- function(lcol){
    if(is.list(lcol)==TRUE){
      ulcol <- unlist(lapply(lcol,paste,collapse=","))
    }
    return(ulcol)
  }
  
  motifs_cytoscape <- as.data.frame(apply(motifs_cytoscape,2,ul),stringsAsFactors = F)
  
  splitmot <- unique(unlist(strsplit(motifs_cytoscape$motif, ",")))
  
  mat <- matrix("0.00",nrow(motifs_cytoscape),length(splitmot))
  colnames(mat) <- splitmot
  
  for (i in 1:nrow(motifs_cytoscape)){
    w <- match(unlist(strsplit(motifs_cytoscape$motif[i],",")), colnames(mat))
    mat[i,w] <- unlist(strsplit(motifs_cytoscape$overlap[i],","))
  }
  
  mat <- cbind(motifs_cytoscape,mat)
  colnames(mat)[1] <- "scans"
  
  return(list(edges = edges, nodes = mat))
}