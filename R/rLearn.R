#' Predict cell type of scATAC-seq dataset by rLearn
#'
#' @param counts A peak-cell count matrix or Seurat object of scATAC-seq data.
#' @param markerfile A list contains the marker region of each cell type.
#'
#' @return A vector of predict label.
#'
#' @import Seurat
#'
#' @export
#'
#'
rLearn = function(counts,
                  markers
){

  marker_count = as.matrix(counts[unique(unlist(markers)),])

  hcc = GetHCC(marker_count, markers, t1=0.5, t2=0.8)

  F = LGC(as.matrix(meta[,9:10]), markers, hcc, weight=T)

  label = Unlabeled(F, markers, t3=0.2)

  result = list()

  result$F = F
  result$label = label

  return(result)

}


#' Read marker from plain text, each row represent a celltype, the first
#' column is the tissue, the second column is celltype name, the latter is marker gene.
#'
#' @param con A connection object or a character string.
#' @param sep The field separator character.
#' @param ... The parameter used to readLines.
#'
#' @return Returns a list of marker genes.
#'
#' @export
#'
ReadMarker = function(con,
                      sep = '\t',
                      ...){

  markerlist = strsplit(readLines(con = con, ...), split = sep)
  names(markerlist) = lapply(markerlist,function(x) x[1])

  markerlist = lapply(markerlist,function(x){
    x = unique(x[-1])
    null_index = which(x == '')
    if(length(null_index) > 0){
      x = x[-null_index]
    }
    return(x)
  })

  return(markerlist)

}

#' Get high confidence cells.
#'
#' @param traindata The scaled data of scRNA-seq.
#' @param markerfile A list contains the marker gene of each cell type.
#' @param t1 The marker gene expression threshold to determine high confidence cells, default 0.5.
#' @param t2 The correlation threshold to determine high confidence cells, default 0.8.
#' @param cor_method Method to calculate correlation, eg. proxyC and base. Default: proxyC.
#'
#' @import outliers proxyC
#'
#' @return A vector of high confidence cells.
#'
#' @export
#'
#'
GetHCC = function(traindata,
                  markerfile,
                  t1,
                  t2,
                  method = 'proxyC'){

  HCC = rep(0, ncol(traindata))

  for(k in 1:length(markerfile)){

    marker_data = as.data.frame(traindata[which(rownames(traindata) %in% markerfile[[k]]),])

    if(nrow(marker_data) == 0){next}

    # threthold = max(colMeans(marker_data)) * t1
    threthold = quantile(colMeans(marker_data),t1)

    filter_cell = names(which(colMeans(marker_data) > threthold))

    if(length(filter_cell) <= 2){next}

    allmarkerdata = traindata[which(rownames(traindata) %in% unlist(markerfile)), filter_cell]

    allmarkermatrix = markerconvert(markerfile)

    allmarkermatrix = allmarkermatrix[rownames(allmarkerdata),]

    if(method == 'proxyC'){

      cor_data = simil(t(as.matrix(allmarkermatrix)),
                       t(as.matrix(allmarkerdata)),
                       method = "cosine",
                       drop0=TRUE)

    }else{

      cor_data = apply(allmarkerdata, 2, function(y){
        apply(allmarkermatrix, 2, function(x) cor(x,y))
      })

    }

    celltype_max = apply(cor_data,2, which.max)

    if(sum(celltype_max == k) == 0){

      next

    }else{

      if(length(markerfile[[k]]) == 1){

        corr_threshold = max(cor_data[k, celltype_max == k]) * t2
        corr_cell = apply(cor_data, 2, function(x) max(x, na.rm = T))
        high_quality_cell = colnames(cor_data)[celltype_max == k &  corr_cell > corr_threshold]

        if(length(high_quality_cell) == 1){

          exp = t(data.frame(high_quality_cell = sapply(markerfile,function(x){mean(traindata[intersect(rownames(traindata), x), high_quality_cell])})))

        }else{

          exp = sapply(markerfile,function(x){

            apply(traindata[intersect(rownames(traindata), x), high_quality_cell], 2, mean)

          })

        }

        dixon_result = apply(exp,1,function(x) dixon.test(x,two.sided = F))
        celltype = sapply(dixon_result,function(x) gsub('^Q\\.','',names(x[[1]])))
        pvalue = sapply(dixon_result,function(x) x$p.value)
        high_quality_cell = high_quality_cell[celltype == names(markerfile)[[k]] & pvalue == min(pvalue)]

      }else{

        corr_threshold = max(cor_data[k, celltype_max == k]) * t2
        corr_cell = apply(cor_data, 2, function(x) max(x, na.rm = T))
        high_quality_cell = colnames(cor_data)[celltype_max == k &  corr_cell > corr_threshold]

      }

    }

    HCC[which(colnames(traindata) %in% high_quality_cell)] = k

  }

  return(HCC)
}

#' marker list to matrix
#'
#' @param markers A list contains the marker gene of each cell type, the list name is cell type.
#'
#' @return A 0-1 matrix that row is marker gene, col is cell type.
#'
#' @export
#'
#'
markerconvert = function(markers){

  genes = sort(unique(unlist(markers)))
  marker_matrix = as.data.frame(matrix(0, length(genes), length(markers)))
  rownames(marker_matrix) = genes
  colnames(marker_matrix) = names(markers)

  for(i in 1:length(markers)){
    marker_matrix[markers[[i]], i] = 1
  }

  return(marker_matrix)

}

#' Get high confidence cells.
#'
#' @param embedding The reduced dimensional data of scRNA-seq.
#' @param markerfile A list contains the marker gene of each cell type.
#' @param HCC A vector of high confidence cells.
#' @param weight A binary value whether to weight high confidence cells.
#'
#' @importFrom parallelDist parallelDist
#'
#' @return A vector of high confidence cells.
#'
#' @export
#'
#'
LGC = function(embedding,
               markerfile,
               HCC,
               weight){

  dm = as.matrix(parallelDist(embedding))

  W = exp(-dm / (2 * 0.2 ^ 2))
  diag(W) = 0

  D = sqrt(matrix(rep(apply(W, 1, sum), nrow(W)),
                  nrow = nrow(W),
                  ncol = ncol(W),
                  byrow = T) * apply(W, 1, sum))
  S = W / D
  S[is.nan(S)] = 0
  S[is.infinite(S)] = 0

  Y = matrix(0,
             nrow = nrow(W),
             ncol = length(markerfile))
  colnames(Y) = 1:(length(markerfile))

  if(weight){

    HCC_weight = log(table(HCC[HCC>0])+1) / max(log(table(HCC[HCC>0])+1))

    for(i in 1:length(HCC)){
      Y[i, which(colnames(Y) == HCC[i])] = 1 * HCC_weight[as.character(HCC[i])]
    }

  }else{

    for(i in 1:length(HCC)){
      Y[i, which(colnames(Y) == HCC[i])] = 1
    }

  }

  I = matrix(0,
             nrow = nrow(S),
             ncol = ncol(S))

  diag(I) = 1

  # F = crossprod(solve(I - 0.99 * S), Y)

  F = SFC:::eigenMapMatMult(S,Y) * 0.99 + (1 - 0.99) * Y
  for(i in 1:20){
    print(i)
    F = SFC:::eigenMapMatMult(S,F) * 0.99 + (1 - 0.99) * Y
  }

  return(F)

}

#' Set unlabeled cells.
#'
#' @param F The predict score matrix of each cells (cells * celltype).
#' @param markerfile A list contains the marker gene of each cell type.
#' @param t3 The threshold to determine the unlabeld cells.
#'
#' @return A predict vector with unlabeled cells.
#'
#' @export
#'
#'
Unlabeled = function(F,
                     markerfile,
                     t3){

  probability = (1/log(F))/(rowSums(1/log(F)))

  max_probability = apply(probability, 1, max)
  max_index = apply(probability, 1, which.max)

  high_probability = c()
  for(i in 1:ncol(probability)){
    high_probability = c(high_probability,quantile(probability[apply(probability,1,which.max) == i, i], 0.9))
    # high_probability = c(high_probability,max(probability[apply(probability,1,which.max) == i, i]) * 0.9)
  }

  F_unlabel = names(markerfile)[max_index]

  F_unlabel[max_probability < high_probability[max_index] * t3] = 'unlabeled'

  return(F_unlabel)

}

#' Evaluate the accurate of hcc.
#'
#' @param hcc A vector of high confidence cells.
#' @param markerfile A list contains the marker gene of each cell type.
#' @param truelabel A vector of true cell type.
#'
#' @return A vector of accuracy of each cell type.
#'
#' @export
#'
#'
hcc_accuracy_estimate = function(hcc,
                                 markerfile,
                                 truelabel){
  acc = c()
  for(i in 1:length(markerfile)){
    acc = c(acc, sum(truelabel[which(hcc==i)] == names(markerfile)[i])/length(which(hcc==i)))
  }
  return(acc)
}

#' load R object and assign new name.
#'
#' @param fileName File path of RData.
#'
#' @return R object.
#'
#' @export
#'
loadRData = function(fileName){

  load(fileName)
  get(ls()[ls() != "fileName"])

}








