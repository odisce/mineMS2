#' @include references.R
#' @include MzDiffFormula.R
#' @include fragPattern.R

#####Labels.R

getMinEdgeLabel <- function(g,v,lab="weight",maxsol=TRUE){
  ev <- incident(g,v,"in")
  vattr <- edge_attr(g,ev,name = lab)
  # vsol <- edge_attr(g,ev,name="solo")
  return(ev[which.min(vattr)])
  
}

getSubformulaLossVertices_afg <- function(g,occs,dags,edge_labels,atoms,label_origin,edge_label="lab"){
  
  formula_labels <- edge_labels$formula
  elabs <- edge_attr(g,"lab")
  uelabs <- unique(elabs)
  idm <- match(elabs,uelabs)
  
  rlf <- sapply(formula_labels[uelabs],MzDiffFormulaFromSingleString,ref=atoms,sep="|")
  rlf_solo <- sapply(rlf,isKnown)
  
  #####Function body.
  V_formula <- matrix(0,nrow=vcount(g),ncol=length(atoms))
  colnames(V_formula) <- colnames(rlf[[1]]@formula)
  
  ###Store the neutral losses corresponding to the edges.
  tempg <- induced_subgraph(g,V(g))
  mdiff <- getMzDiff(tempg,occs = occs,dags,edge_labels$mz)
  solv <- rlf_solo[match(edge_attr(tempg,"lab"),uelabs)]
  vmloss <- 1/(1+apply(mdiff,2,mean))+ifelse(solv,0,2000)
  
  
  tempg <- set_edge_attr(tempg,"weight",value = vmloss)
  c_vertex <- 1
  
  ####We get a V_sequence reordered by the label origin.
  e_out <- incident(tempg,1)
  v_seq <- head_of(tempg,e_out)
  e_attr <- edge_attr(graph = tempg,name="lab",index = e_out)

  for(iv in seq_along(v_seq)){
    v <- v_seq[iv]
    
    ####
    e_inc <- getMinEdgeLabel(tempg,v)
    if(length(e_inc)==0){
      next
    }else{
      v_source <- tail_of(tempg,e_inc)
      v_source_pos <- match(as.numeric(v_source),as.numeric(V(tempg)))
      e_lab <- edge_attr(graph=tempg,index = e_inc,name = "lab")
      pfor <- which(uelabs==e_lab)
      e_for <- rlf[[pfor]]
      
      ####If there is a single formula we update the information:
      if(isKnown(e_for)){
        ###We update the possible formula of the loss
        V_formula[iv+1,] <- V_formula[v_source_pos,]+e_for@formula[1,,drop=FALSE]
      }else{
        ###nothing change
        V_formula[iv+1,] <- V_formula[v_source_pos,]
      }
    }
  }
  ###We compute the corresponding edge sequence labels.
  return(V_formula)
}

chooseVerticesLosses <- function(lf,oformula,dags,mzv,vrdbe,atoms,subformula = NULL){
  ###Solo
  if(length(lf) == 1){
    if(length(oformula)==0){
      return(c(1,1,0))
    }
    vm <- sapply(oformula,isSubformula,x=lf,simplify=TRUE)
    return(c(1,1,as.numeric(all(vm))))
  }
  
  bloss <- NULL
  bpartial <- 0
  bsolo <- 0
  pok <- 1:length(lf)
  if((length(oformula) == 0) || ((length(oformula) == 1) && is.na(oformula))) {
    bloss <- 1:length(lf)
  } else {
    ####If a possible subformula is furnished, the possible formula of the vertices are restricted.
    if(!is.null(subformula)){
      vm_sub <- isSubformula(as.integer(subformula),lf)
      pok <- which(vm_sub)
      if(length(pok)>0){
        lf <- lf[pok,]
      }
      if(length(pok)==0){
        pok <- 1:length(lf)
        return(c(UNKNOWN_FORMULA_INDICATOR,0,0))
      }
    }
    
    if((length(lf) > 1) && (length(oformula) > 0)) {
      ###Majority vote.
      vm <- sapply(oformula,isSubformula,x=lf,type="none",simplify=TRUE)
      bf <- apply(vm,1,sum)
      bloss <- which(bf == max(bf))
    } else {
      bloss <- 1
    }
  }
  
  # We keep the one the closest to the furnished mass.
  if(length(bloss) > 1) {
    m_atoms <- getAtomMass(atoms)
    vmz <- lf@formula[bloss,] %*% m_atoms
    rdbe <- calcRDBE_raw(lf@formula[bloss,],vrdbe)
    
    # Now we split the value of t
    bloss <- bloss[which.min(abs(vmz-mzv)*10^6/mzv+abs(rdbe))]
    if ((length(lf) > 1) && (length(oformula) > 0)) {
      vm <- sapply(oformula,isSubformula,x=lf,type="none",simplify=TRUE)
      bf <- apply(vm,1,sum)
      bpartial <- as.numeric(as.numeric(bf[bloss]) == nrow(lf@formula))
    }
  } else {
    bsolo <- 1
  }

  # We return the subformula with the maximum size
  resv <- c(pok[bloss], bsolo, bpartial)
  # if(abs((mzv-161.0362)<0.002)) browser()
  
  return(resv)
  
}

###allf is supposed to be a set of MzDiffFormula corresponding to the full set of edge labels.
annotateVertices <- function(fp,vlab,dags,elabs,atoms,allf = NULL,massdiff=NULL, edge_label="lab",oformula=NULL){
  g <- mm2Graph(fp)
  occs <- mm2Occurrences(fp)

  if(is.null(massdiff)){
    massdiff <- getMzDiff(g,occs,dags,elabs$mz)
    massdiff <- apply(massdiff,2,mean)
  }

  subformula <- getSubformulaLossVertices_afg(g,occs,dags,elabs,atoms,c(0,vlab))
  
  if(is.null(allf)){
    allf <- sapply(elabs$formula[vlab],function(x,atoms){
      MzDiffFormulaFromSingleString(x,ref = atoms,sep = "|")
    },atoms=atoms,simplify=TRUE)
  }

  
  vrdbe <- vecRDBE(allf[[1]]@formula)
  mzv <- massdiff[seq_along(vlab)]

  # if(is.null(oformula)){
  #   oformula <- sapply(getFormula(m2l)[occs[,1]],MzDiffFormulaFromSingleString,ref=atoms,sep="|")
  # }
  
  ####Each col correspond to a vertices.
  res <- sapply(seq_along(vlab),function(x,l_origin,voformula,vatoms,subf,vdags,lformula,vmz, vrdbe){
    
    ####We build the formula corresponding to the loss
    lf <- lformula[[x]]
    if(length(lf)==0){
      return(list(c(1,1),lf))
    }
    vrdbe <- vecRDBE(lf@formula)
    rtemp <- chooseVerticesLosses(lf,oformula=voformula,dags,vmz[x], vrdbe=vrdbe,atoms=vatoms,subformula = subf[x+1,])

    return(list(rtemp[2:3],lf[rtemp[1],,drop=FALSE]))
  },vdags=dags,voformula=oformula,l_origin=vlab,vatoms = atoms,subf = subformula,
  lformula=allf,vmz = mzv, vrdbe=vrdbe)
  return(res)
}



#Return the subformula coherent with the most formula.
makeEdgeLabelWithoutVertices <- function(lf,mzv,oformula,atoms,vrdbe){
  
  ###Solo
  if(length(lf)==1){
    vm <- sapply(oformula,isSubformula,x=lf,simplify=TRUE)
    return(c(1,1,as.numeric(all(vm))))
  }
  if(length(lf)==0){
    return(c(NA_real_,1,1))
  }
  
  bloss <- NULL
  bpartial <- 0
  bsolo <- 0
  bf <- numeric(0)
  if((length(oformula)==0) || ((length(oformula)==1) && is.na(oformula))) {
    bloss <- 1:length(lf)
  }else{
    ###Majority vote.
    vm <- sapply(oformula,isSubformula,x=lf,type="none",simplify=TRUE)
    
    bf <- apply(vm,1,sum)
    bloss <- which(bf==max(bf))
  }
  
  ###We keep the one the closest to the furnished mass.
  if(length(bloss)>1){
    m_atoms <- getAtomMass(atoms)
    vmz <- lf@formula[bloss,] %*% m_atoms
    rdbe <- calcRDBE_raw(lf@formula[bloss,],vrdbe)
    ###mono atom.
    penalisation <- apply(lf@formula[bloss,],1,function(x){
      vs <- sort(x[c("C","H","N","O")],decreasing = TRUE)
      (vs[1]-vs[2])/vs[1]
    })
    scorev <- abs(vmz-mzv)*10^6/mzv+abs(rdbe)/10+penalisation
    bloss <- bloss[which.min(scorev)]
    if(length(bf)==0){
      bpartial <- TRUE
    }else{
      bpartial <- as.numeric(as.numeric(bf[bloss])==nrow(lf@formula))
    }
  }else{
    bsolo <- 1
  }
  ###We return the subformula with the maximum size
  c(bloss,bsolo,bpartial)
}

makeEdgeLabel <-
  function(fp,
           vf,
           allf,
           vrdbe,
           atoms,
           oformula,
           vmass) {
    
    g <- mm2Graph(fp)
    e <- E(g)
    
    ####We calculate the putative function.
    mat_dist <-
      matrix(c(tail_of(g, e), head_of(g, e)), byrow = FALSE, ncol = 2)
    coevec <- rep(TRUE,length=ecount(g))
    multipleannot <- rep(FALSE,length=ecount(g))
    
    ###Formula annotation of the edges.
    lf <- vector(mode="list",length=ecount(g))
    
    p_non_root <- which(mat_dist[, 1] != 1)
    p_root <- which(mat_dist[,1]==1)
    
    #####REORDERING TO TAKE INTO ACCOUNT THE WRONG FORMULA
    pm <- match(mat_dist[p_root,2],as.numeric(V(g)))-1
    lf[p_root[pm]] <- vf
    
    vannot <- rep(NA_character_,ecount(g))
    if (length(p_non_root) > 0) {
      
      ###We compute the difference of formula.
      vdist <-
        apply(mat_dist[p_non_root, , drop = FALSE], 1, function(x, vf) {
          fa <- vf[[x[2]-1]]
          fb <- vf[[x[1]-1]]
          if((length(fa)==0)|(length(fb)==0)){
            return(rep(NA_integer_,ncol(fa@formula))) 
          }
          fa@formula - fb@formula
        },vf=vf)
      
      
      vdist <- unname(split(vdist, rep(1:ncol(vdist), each = nrow(vdist))))
      
      ###Mass with NA labels.
      to_compute <-sapply(vdist,function(x) !is.na(x[1]))
      
      ###Now we check if this formula is in the corresponding labels.
      rv <- vector(mode="list",length=length(vdist))
      for(i in seq_along(rv)){
        rv[[i]] <- numeric(0)
      }
      
      if(length(to_compute)>0){
        rv[to_compute] <- mapply(
          vdist[to_compute],
          allf[p_non_root[to_compute]],
          FUN = function(x, y) {
            idxFormula(x, y)
          },
          SIMPLIFY = FALSE
        )
      }
      
      
      ####Coherent and non coherent edge
      coherent_edge <- sapply(rv, length) > 0
      non_coherent_edge <- !coherent_edge
      high_mz <- non_coherent_edge&(!to_compute)
      coherent_edge <- which(coherent_edge)
      non_coherent_edge <- which(non_coherent_edge)
      
      if (length(coherent_edge) > 0) {
        vannot[p_non_root[coherent_edge]] <- rv[coherent_edge]
        ###We process the rest of the formula.
        for(cp in coherent_edge){
          lf[[p_non_root[cp]]] <- allf[[p_non_root[cp]]][rv[[cp]],drop=FALSE]
        }
        coevec[p_non_root[coherent_edge]] <- TRUE
      }
      
      if (length(non_coherent_edge) > 0) {
        ###mzv is the mass of the corresponding loss.
        #oformula correspond to the et of molecular formula of the occurrences.
        #atoms id the names of the atoms.
        #vrdbe represent the rdbe header.
        tempe <- mapply(
          allf[p_non_root[non_coherent_edge]],
          vmass[p_non_root[non_coherent_edge]],
          FUN = makeEdgeLabelWithoutVertices,
          MoreArgs = list(
            'vrdbe' = vrdbe,
            'atoms' = atoms,
            'oformula' = oformula
          )
        )
        pna <- is.na(tempe[1,])
        pnna <- !pna
        vannot[p_non_root[non_coherent_edge[pna]]] <- HIGH_MASS_INDICATOR
        vannot[p_non_root[non_coherent_edge[pnna]]] <- tempe[1,pnna]
        coevec[p_non_root[non_coherent_edge]] <- FALSE
        coevec[p_non_root[high_mz]] <- TRUE
        multipleannot[p_non_root[non_coherent_edge[!as.logical(tempe[2,])]]] <- TRUE
        for(cp in which(pnna)){
          lf[[p_non_root[non_coherent_edge[cp]]]] <- allf[[p_non_root[non_coherent_edge[cp]]]][tempe[1,cp],drop=FALSE]
        }
        if(sum(pna)>0)
          for(cp in which(pna)){
            lf[[p_non_root[non_coherent_edge[cp]]]] <- MzDiffFormula(ref=atoms)
          }
       }
    }
    return(list(formula=lf,annot=as.integer(vannot),root=ifelse(mat_dist[,1]==1,TRUE,FALSE),
                coherent = coevec, multiple = multipleannot))
  }

#' All the exact m/z differences for a given fragPattern object
#'
#' Return all the exact m/z differences for a given fragPattern object.
#' Should not be used by the user.
#'
#' @param g fragPattern object
#' @param occs list of occurrences in the fragPattern
#' @param mgs list of fragmentation graphs (DAGs) contained in the fragPattern
#' @param loss_mz list of all discreized the m/z differences
#'
#' @export
getMzDiff <- function(g, occs, mgs, loss_mz){
  occs_gid <- occs[, 1]
  occs_pos <- occs[, 2]
  res_diff <- matrix(ncol=ecount(g),nrow=nrow(occs))
  for (i in 1:length(occs_gid)) {
    gid <- occs_gid[i]
    pos <- occs_pos[i]
    map <- get_mapping(mgs[[gid]], g, loss_mz, root = pos)
    mzs <- vertex_attr(mgs[[gid]], "mz")
    ids <- V(mgs[[gid]])
    ###Matched peaks.
    matched_peaks_idx <- match(map[2,], ids)
    
    ###
    V_head <- head_of(g,es = E(g))
    V_tail <- tail_of(g,es = E(g))
    posh <- match(V_head,map[1,])
    post <- match(V_tail,map[1,])
    
    diffv <- mzs[matched_peaks_idx[post]]-
      mzs[matched_peaks_idx[posh]]
    res_diff[i,] <- diffv
  }
  
  
  return(res_diff)
}

#' Calculate labels (m/z differences) for vertices and edges of a pattern
#' 
#' For a specific pattern, this function calculates the labels (m/z differences) for vertices and edges.
#' Labels for the vertices correspond to m/z differences from the root of the pattern (precursor parent).
#' 
#' This function should not be used by the user.
#' @param fp id of the pattern
#' @param atoms vector of atom names to consider (for example, CHNOPS)
#' @param dags vector of the fragmentation graphs (DAGs) contained in the pattern
#' @param elabs dataframe of all the discretized edge labels (m/z differences) in the dataset
#' @param oformula vector of precursor formula for each spectra in the pattern, if given
#' @param edge_label id of the field containing the m/z differences in elabs
#' 
#' @export
annotateAFG <- function(fp,atoms,dags,elabs,oformula,edge_label="lab"){
  ###We first build a set of LossGraph corresponding to all the edge labels.
  g <- mm2Graph(fp)
  occs <- mm2Occurrences(fp)
  atoms <- atoms
  
  ####Getting the vertices and the edges to pass, in the right order
  all_edges <- E(g)
  p_root <- tail_of(g,all_edges)==1
  e_origin <- all_edges[p_root]
  v_seq <- head_of(g,es = e_origin)
  
  pmm <- match(as.numeric(v_seq),as.numeric(V(g)))-1
  
  ###LABEL ORIGIN IS IN THE SAME ORDER THAN V, DO NOT MOVE.
  label_origin <- edge_attr(g,name=edge_label,index = e_origin)[pmm]
  
  ###Building all the formula.
  all_edges_lab <- edge_attr(g,name = "lab")
  allf <- sapply(elabs$formula[all_edges_lab],function(x,atoms){
    MzDiffFormulaFromSingleString(x,ref = atoms,sep = "|")
  },atoms=atoms,simplify=FALSE)
  
  massdiff <- getMzDiff(g,mm2Occurrences(fp),dags,elabs$mz)
  massdiff <- apply(massdiff,2,mean)

  ## mean values of mz diff on the full dataset
  massdiff_ref <- elabs[edge_attr(g, name="lab"),"mz"]
  
  ###If there is an NA we replace the corresponding mz value by the lab value.
  pna <- is.na(massdiff)
  if(sum(pna)>0){
    ###We get the corresponding value in the lab table?
    massdiff[pna] <- elabs$mz[all_edges[pna]]
  }

  vannot <- annotateVertices(fp,label_origin,dags,elabs,atoms,allf = allf[pmm],
                             edge_label="lab",oformula = oformula,massdiff = massdiff[pmm])
  formula_vertices <- vector(mode="list",length=ncol(vannot))
  
  df_vertices <- do.call(rbind,vannot[1,])

  for(i in seq_along(formula_vertices)){
    formula_vertices[i] <- vannot[2,i]
  }
  
  ####FORMULA VERTICES IS IN THE RIGHT ORDER DON'T TOUCH IT.
  elab <- makeEdgeLabel(fp,
                        vf = formula_vertices,
                        allf = allf,
                        vrdbe = vecRDBE(allf[[1]]@formula),
                        atoms = atoms,
                        oformula = oformula,
                        vmass = massdiff) 
  #return(list(vertices = list(formula=formula_vertices,info=df_vertices),edges=elab,
  #            mz=list(v=massdiff[p_root][pmm],e=massdiff)))
  return(list(vertices = list(formula=formula_vertices,info=df_vertices),edges=elab,
              mz=list(v=massdiff_ref[p_root][pmm],e=massdiff_ref)))
}

makeLabelPattern <- function(p,atoms,dags,elabs,oformula,mzdigit=4, print_formula=TRUE){
  annot <- annotateAFG(p,atoms,dags,elabs,oformula = oformula)
  g <- mm2Graph(p)
  
  labvert <- character(vcount(g)-1)
  labedge <- character(ecount(g))
  
  ####We build the vertices label.
  vert <- annot$vertices
  labvert <- paste("PP -",sprintf(paste("%0.",mzdigit,"f",sep=""),
                                 annot$mz$v),sep=" ")
  
  labvert_formula <- sapply(vert$formula,function(x,def){
    if(length(x)==0) return(def)
    as.character(x)
  },def=HIGH_MASS_LEGEND)
  
  ####If the formula is not solo or partial we add parenthesis
  vpv <- apply(!vert$info[,1,drop=FALSE],1,any) ## vpv : logical vector : if TRUE : several formula

  
  labvert_formula[vpv] <- paste("(",labvert_formula[vpv],")",sep="")  ## add parenthesis
  
  if(print_formula) labvert <- c("PP",paste(labvert,labvert_formula,sep="\n")) ## to add formula to label edges
  else labvert <- c("PP",labvert)
  
  
  #### We build the edges label.
  edg <- annot$edges
  
  labedge <- sapply(edg$formula,function(x,def){
    if(length(x)==0) return(def)
    as.character(x)
  },def=HIGH_MASS_LEGEND)
  
  vp <- ((!edg$coherent)&(labedge!=HIGH_MASS_LEGEND))
  
  vp[1:(vcount(g)-1)] <- vpv
  pmz <-  labedge==HIGH_MASS_LEGEND
  labedge[pmz] <- sprintf(paste("%0.",mzdigit,"f",sep=""),annot$mz$e[pmz])

  ## if several formula, put parenthesis or mz
  vp <- vp|((edg$multiple)&(!pmz))
  
  if(print_formula)
  {
    ## only one formula
    labedge[!vp&(!pmz)] <- paste(sprintf(paste("%0.",mzdigit,"f",sep=""),annot$mz$e[!vp&(!pmz)]), labedge[!vp&(!pmz)], sep="\n") ## mz + formula 

    #several formula
    labedge[vp] <- paste(sprintf(paste("%0.",mzdigit,"f",sep=""),annot$mz$e[vp]), paste("(",labedge[vp],")",sep=""), sep="\n") ## mz + parenthesis 
  }
  else{
  ## mz for all edges
    labedge <- sprintf(paste("%0.",mzdigit,"f",sep=""),annot$mz$e)
  }

  ####For all the high mass label we put the mz diff over the formula.
  return(list(vertices=labvert,edges=labedge))
}
