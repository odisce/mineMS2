#' @include references.R
#' @include MzDiffFormula.R

MULTIPLE_FORMULA <- "Formula>2"
HIGH_MASS <- "high mz"

###Function used to construct the edge labels dataset.
make_label_loss <- function(lab_edges){
  
  formula <- lab_edges$formula
  
  str_formula <- str_split(formula,pattern = fixed("|"))
  
  str_formula <- sapply(str_formula,function(x){
    
    ###If no label have been found.
    if(anyNA(x)) return(HIGH_MASS)
    if(length(x)==1) return(x)
    if(length(x)>1) return(MULTIPLE_FORMULA)
  })
  
  full_labs <-  paste("- ",sprintf("%0.3f",lab_edges$mz),"\n",str_formula,sep="")
  
  
  labs <- ifelse(((str_formula==HIGH_MASS)|(str_formula==MULTIPLE_FORMULA)),
                 sprintf("%0.3f",lab_edges$mz),str_formula)
  return(data.frame(labs=labs,full_labs=full_labs))
}

make_label_frag <- function(lab_frags){
  
  formula <- lab_frags$formula
  
  str_formula <- str_split(formula,pattern = fixed("|"))
  
  str_formula <- sapply(str_formula,function(x){
    
    ###If no label have been found.
    if(anyNA(x)) return(HIGH_MASS)
    if(length(x)==1) return(x)
    if(length(x)>1) return(MULTIPLE_FORMULA)
  })
  
  full_labs <-  paste(sprintf("%0.3f",lab_frags$mz),"\n",str_formula,sep="")
  
  
  labs <- ifelse(((str_formula==HIGH_MASS)|(str_formula==MULTIPLE_FORMULA)),
                 sprintf("%0.3f",lab_frags$mz),str_formula)
  return(data.frame(labs=labs,full_labs=full_labs))
}



reorderAtom <- function(atoms){
  
  lf<- MzDiffFormula(ref=names(atoms))
  return(atoms[match(colnames(lf@formula),names(atoms))])
}

#' Discretize m/z differences.
#'
#' It also builds fragmentation graphs from MS/MS spectra using discretized m/z differences as edges.
#'
#' @param m2l The ms2 lib object to be discretized.
#' @param ppm the maximum authorized deviation in ppm (parts per million) [default: 15]
#' @param dmz The maximum authorized deviation in Da [default: 0.007]
#' @param count The minimum number of spectra in which the label needs to be found. Shall be greater or equal to 2 [default: 2]
#' @param limMzFormula An interval giving the range in which the formula will be calculated: peaks with a m/z lower than the lower term will be ignored while peak with a mass higher than the higher term won't have formula check [default: c(14.5, 200)]
#' @param maxFrags The maximum number of fragment allowed on the spectra [default: 15]
#' @param maxOverlap The degree of overlap allowed between bins. Bins with a higher overlap are fused [default: 0.05]
#' @param strictMatching Shall the formula matching be strict, or approximated [default: TRUE]
#' @param precPpm The ppm tolerance used to match the precursor mz to the fragments [default: 20]
#' @param precDmz The minimum tolerance used to match the precursor mz to the fragments [default: 0.02]
#' @param atoms A list of the maximum number of atoms authorized during the formula generation process. If set to NULL [default], the default list includes C:max(limMzFormula)/12, H:50, N:10, O:15, Cl:2, S:2, P:2 (if heteroAtoms is set to its default, TRUE), or C:max(limMzFormula)/12, H:50, N:10, O:15 otherwise
#' @param heteroAtoms Should heteroatoms P, Cl, S be included as possible atoms (ignored if a custom list of atoms is provided) [default: TRUE]
#' @param mzDigits The number of decimals to display on the graphs [default: 4]
#'
#' @return An ms2Lib object with filled fields and constructed dags.
#' @export
#'  
#' @rdname discretizeMzDifferences
#' 
#' @examples
#' # Loading the data
#' data(m2l)
#' 
#' # Constructing edge labels
#' m2l <- discretizeMzDifferences(m2l)
#' 
#' # Setting the maximum number of atoms to specific values
#' m2l <- discretizeMzDifferences(m2l, 
#' atoms = list("C"=50,"H"=100,"N"=6,"O"=6,"S"=2,"Cl"=1,"P"=2))
setMethod("discretizeMzDifferences", "ms2Lib", function(m2l,
                                                        ppm = 15,
                                                        dmz = 0.007,
                                                        count = 2,
                                                        limMzFormula = c(14.5, 200),
                                                        maxFrags = 15,
                                                        maxOverlap = 0.05,
                                                        strictMatching = TRUE,
                                                        precPpm = 20,
                                                        precDmz = 0.02,
                                                        atoms = NULL,
                                                        heteroAtoms = TRUE,
                                                        mzDigits = 4) {
  
  message("Discretization of the m/z differences...")
  
  ###Parameters checking
  if(maxFrags>20){
    stop("mineMS2 is not made to process more than 20 fragments by spectra, please lower
			 the maxFrags parameters")
  }
  if(maxOverlap>0.15){
    warning("maxOverlap is superior to 0.15, this can lead to incoherent labels.")
  }
  
  if((!is.numeric(limMzFormula))|(length(limMzFormula)!=2)){
    stop("limMzFormula should be a range with 2 floats")
  }
  
  if((limMzFormula[1]<=10|limMzFormula[2]>250)){
    limMzFormula <- c(max(limMzFormula[1],10),min(max(limMzFormula[2],250)))
    warning("limMzFormula is to wide, to avoid computation issue it have been set to ",
            limMzFormula[1],"-",limMzFormula[2])
  }
  
  if(count==1){
    warning("count set to a 2.")
    count <- 2
  }
  freq <- count/length(m2l@spectra)
  if((freq<0) | (freq>1)){
    stop("Wrong count value ",count)
  }
  
  if(is.null(atoms)){
    if(heteroAtoms){
      #atoms <- list("C"=max(limMzFormula)%/%12,"H"=50,"N"=6,"O"=6,"S"=2,"Cl"=1,"P"=2)
      atoms <- list("C"=max(limMzFormula)%/%12,"H"=50,"N"=10,"O"=15,"Cl"=2,"S"=2,"P"=2)
    }else{
      atoms <- list("C"= max(limMzFormula)%/%12,"H"=50,"N"=10,"O"=15)
    }
  }
  
  m2l@atoms <- reorderAtom(atoms)
  res_list <- .discretizeMzDifferences(mm2Spectra(m2l),
                                          ppm = ppm, dmz = dmz,
                                          freq = freq, mzdigits = mzDigits,
                                          limFormula = limMzFormula, maxFrag = maxFrags,
                                          max_overlap = maxOverlap, strictMatching = strictMatching,
                                          prec.ppm = precPpm,prec.dmz = precDmz, atoms = m2l@atoms)
  
  
  ###At this step we add the formula if necessary.
  
  ###We find the position of the last elementsapply()
  plast <- which(res_list$elems$mz<limMzFormula[2])
  plast <- plast[length(plast)]
  ref <- names(atoms)
  
  message("Formula extensions")
  time1 <- Sys.time()
  to_correct <- find_combinations_ranges(res_list$elems[1:(plast+1),"mzmin"],res_list$elems[1:(plast+1),"mzmax"],limMzFormula[2])
  print(Sys.time() - time1)

  time2 <- Sys.time()
  ###We merge the formula when necessary.
  allF <- sapply(res_list$elems$formula[1:plast],function(x,atoms){
    orderByRDBE(MzDiffFormula(str_split(x,fixed("|"),simplify=TRUE)[1,],ref=atoms))
  },atoms=names(m2l@atoms))
  
  for(i in seq_along(allF)){
    
    pf1 <- to_correct[[i]]$f1+1
    pf2 <- to_correct[[i]]$f2+1
    
    for(j in seq_along(pf1)){
      
      genf <- combineMzDiffFormula(allF[[pf1[j]]],allF[[pf2[j]]])
      allF[[i]] <- addFormula(genf,allF[[i]])
      
    }
    allF[[i]] <- orderByRDBE(allF[[i]])
  }
  
  ####Now we rebuild the edge labels.
  res_list$elems$formula[1:plast] <- sapply(allF,function(x){
    paste(as.character(x),collapse = "|")
  })
  
  print((Sys.time() - time2))
  ###Add the fusing part of the edge labels.
  message("Fusing the edge labels")
  time1 <- Sys.time()
  ###Simplifying the DAG if necessary.
  change <- TRUE
  niter <- 1
  while(change&niter<4 && nrow(res_list$elems) > 1){
    res_list <- fuseElem(elems =res_list$elems ,dags = res_list$dags,thresh=count,atoms = names(atoms))
    change <- res_list$change
    niter <- niter+1
  }
  
  ###Constructing the edge labels
  templabs <- make_label_loss(res_list$elems)
  res_list$elem$full_labels <- templabs$full_labs
  res_list$elem$labs <- templabs$labs
  
  mm2EdgesLabels(m2l) <- res_list$elems
  mm2Dags(m2l) <- res_list$dags
  print((Sys.time() - time1))
  message("Discretization completed: ", nrow(res_list$elem)," common m/z differences found.")
  m2l
})


removeBinsOverlap <- function(mzmin,mzmax,margin=0.00001){
  tmzmin <- mzmin
  tmzmax <- mzmax
  poverlap <- which(mzmax[1:(length(mzmax)-1)]>mzmin[2:(length(mzmax))])
  if(length(poverlap)>0){
    dev <- (tmzmin[poverlap]+tmzmax[poverlap])/2
    tmzmin[poverlap+1] <- dev-margin
    tmzmax[poverlap] <- dev+margin
  }
  return(list(min=tmzmin,max=tmzmax))
}

.discretizeMzDifferences <- function(list_spec,
                                         ppm = 7, dmz = 0.002,
                                         freq = 0.05, limFormula=c(14.5,200),
                                         mzdigits=3, maxFrag = 25,
                                         max_overlap=0.05,strictMatching=TRUE,
                                         prec.ppm=20,prec.dmz=0.02,
                                         atoms = list("C"= 1000%/%12,"H"=50,"N"=6,"O"=6),...) {
  list_spec <- lapply(list_spec,getExtendedSpec,maxNum=maxFrag,ppm=prec.ppm,
                      dmz=prec.dmz,relInt = TRUE)
  list_masses <- lapply(list_spec,'[[',i="mz")
  list_int <- lapply(list_spec,'[[',i="rel_int")
  
  mean_mzf <- sum(limFormula)/2
  tol <- diff(limFormula)/2
  l_atoms <- atoms
  
  ndiff <- length(list_masses)
  thresh <- checkFracParam(freq, ndiff)
  list_matrix <- lapply(list_masses,generateMatDiff)
  pcord <- lapply(list_matrix, '[[', i = "idx")
  vecsample <- rep(seq_along(pcord), times = sapply(pcord, nrow))
  vecsize <- sapply(pcord, nrow)
  pcord <- do.call('rbind', pcord)
  vmz <- unlist(sapply(list_matrix, '[[', i = "mz"))
  allEdges<- NULL
  valinter <- c(0, cumsum(vecsize))
  treatInter <- function(x) {
    pvec <- .bincode(x, valinter, right = TRUE)
    if (length(pvec) == 1) {
      tr <- matrix(c(pvec, pcord[x, ]), nrow = 1)
      colnames(tr) <- c("pvec", "row", "col")
      return(tr)
    }
    cbind(pvec, pcord[x, ])
  }
  
  ###We start by discretizing all the values
  resdisc <-
    discretizeSequenceByClasses(
      ppm = ppm,
      vmz = vmz,
      vsample = vecsample,
      dmz = dmz,
      ndiff = length(list_masses),
      frac = freq,extended=TRUE
    )
  
  ###Now we create all the standard deviation of the masses.
  masses <- resdisc$discrete_elem[,1]
  diffv <- resdisc$discrete_elem[,3]-resdisc$discrete_elem[,2]
  
  ###We define the max larg
  qv <- resdisc$discrete_elem[,4]
  
  ###The standard deviation is created
  sds <- qv
  
  fac_sig <- 2
  ###Merging overlapping values.
  message("m/z merging")
  
  merged_masses <- gaussianMerging(masses,sds^2,alpha=max_overlap,fac_sig = fac_sig)
  high_mz_idx <- which(merged_masses$mu>limFormula[2])
  ###Now we check the values which may originates from a formula.
  # browser()
  
  ###Generation of all the neutral formula
  message("Formula generation")
  
  ###How ot handle heteroatoms, maximum number of different heteroatoms.
  time1 <- Sys.time()
  allFormula <- lossesFormulaGeneration(limFormula,atoms=l_atoms,...)
  
  
  ##Associating formula with losses.
  f_masses <- allFormula$masses

  ###Now we generate the bmin to fac_sig
  tval <- f_masses*ppm/10^6
  sd_f <- ifelse(tval>dmz,tval,dmz)
  
  maxm <- f_masses+fac_sig*sd_f
  minm <- f_masses-fac_sig*sd_f

  merginf <- merged_masses$mu-fac_sig*sqrt(merged_masses$sig)
  mergmax <- merged_masses$mu+fac_sig*sqrt(merged_masses$sig)
  
  ###We remove the remaining overlaps.
  tempmz <- removeBinsOverlap(merginf,mergmax,margin=0)
  merginf <- tempmz$min
  mergmax <- tempmz$max

  ###Now we need to check the overlap between the interval.
  ##Function to parse the Cpp function.
  #print(minm)
  #print(maxm)
  #print(merginf)
  #print(mergmax)
  matched_inter <- checkInter(minm,maxm,merginf,mergmax)
  vec_inter <- do.call("cbind",matched_inter)
  
  atm_names <- colnames(allFormula$formula)
  ###Now we create a table centralizing all the values.
  str_formula <- NULL
  freq_f <- NULL
  non_freq_f <- NULL
  
  xp1 <- 1
  xp2 <- 2
  if(strictMatching){
    xp1 <- 3
    xp2 <- 4
  }
  print((Sys.time() - time1))

  str_formula <- apply(vec_inter,1,function(x){
    if(is.na(x[xp1])) return(NA_character_)
    apply(allFormula$formula[x[xp1]:x[xp2],,drop=FALSE],1,formulaToString,vnames=atm_names)
  })
  
  ###Final strings
  str_formula <- sapply(str_formula,function(x){
    if((length(x)==1) && is.na(x)) return(NA_character_)
    paste(x,collapse = "|")
  })
  
  ####We check if the losses are nitrogen or carbon only.
  carb_only <- str_detect(str_formula,"^C[0-9]{0,10}$")
  nitrogen_only <- str_detect(str_formula,"^N[0-9]{0,10}$")
  CH_only <- str_detect(str_formula,"^C[0-9]{0,10}H[0-9]{0,10}$")
  
  ###we filter out intervals without formula.
  to_keep <- which(!is.na(str_formula))
  
  ###Now we build a recapitulative table, keeping only the mean mz,min mz, formula and label.
  recap_tab <- data.frame(lab=1:(length(to_keep)+length(high_mz_idx)),mz=merged_masses$mu[c(to_keep,high_mz_idx)],
                          mzmin=merginf[c(to_keep,high_mz_idx)],
                          mzmax=mergmax[c(to_keep,high_mz_idx)],
                          sig=merged_masses$sig[c(to_keep,high_mz_idx)],
                          fused = merged_masses$fused[c(to_keep,high_mz_idx)],
                          count=numeric((length(to_keep)+length(high_mz_idx))),
                          formula=c(str_formula[to_keep],rep(NA_character_,length(high_mz_idx))),
                          carb_only = c(carb_only[to_keep],rep(FALSE,length(high_mz_idx))),
                          nitrogen_only = c(nitrogen_only[to_keep],rep(FALSE,length(high_mz_idx))),
                          ch_only = c(CH_only[to_keep],rep(FALSE,length(high_mz_idx))),
                          stringsAsFactors = FALSE)
  
  ###Graph creations
  toReturn <- vector(mode = "list", length = length(list_masses))
  
  for (i in seq_along(list_masses)){
    toReturn[[i]] <- matrix(
      NA_real_,
      nrow = length(list_masses[[i]]),
      ncol = length(list_masses[[i]])
    )
    colnames(toReturn[[i]]) <- sprintf(paste('%0.',mzdigits,'f',sep=""),list_masses[[i]])
    ###We first need to reorder the mass diff
    om_diff <- order(list_matrix[[i]]$mz,decreasing = FALSE)
    list_matrix[[i]]$mz <- list_matrix[[i]]$mz[om_diff]
    list_matrix[[i]]$idx <- list_matrix[[i]]$idx[om_diff,,drop=FALSE]
    
    matched <- disjointBins(list_matrix[[i]]$mz, recap_tab$mzmin, recap_tab$mzmax, recap_tab$mz)

    ###We keep only the matched edges.
    pok <- which(!is.na(matched))
    
    ###We keep only the labeled edges.
    list_matrix[[i]]$mz <- list_matrix[[i]]$mz[pok]
    list_matrix[[i]]$idx <- list_matrix[[i]]$idx[pok,,drop=FALSE]
    toReturn[[i]][list_matrix[[i]]$idx]<-matched[pok]
    
  
    ###We create the graph.
    idxnodes <- match(list_masses[[i]],sort(list_masses[[i]],decreasing = FALSE))-1
    dag <- make_empty_graph(n = nrow(toReturn[[i]]),directed=TRUE)
    dag <- set_vertex_attr(dag,name="mz",value = list_masses[[i]])
    dag <- set_vertex_attr(dag,name="rel_int",value = list_int[[i]])
    dag <- set_vertex_attr(dag,name="label",value = as.integer(idxnodes))
    dag <- set_vertex_attr(dag,name="prec",value=as.logical(list_spec[[i]]$prec))
    pf <- which(!is.na(toReturn[[i]]),arr.ind = TRUE)
    to_remove <- numeric(100)
    
    nrm <- 1
    if(nrow(pf)!=0){
      ##We remove the parent with the same descendant.
      lab_entry <- paste(pf[,2],toReturn[[i]][pf],sep="_")
      
      pdup <- table(lab_entry)
      
      to_check <- which(as.numeric(pdup)>=2)
      if(length(to_check)!=0){
        for(j in to_check){
          if(nrm>length(to_remove)) ro_remove <- c(to_remove,numeric(length(to_remove)))
          ###We get the concerned edges
          posf <- which(lab_entry==names(pdup)[j])
          prob_edges <- pf[posf,]
          
          ###Now we calculate the masses.
          dmz <- list_masses[[i]][prob_edges[,1]]-list_masses[[i]][prob_edges[,2]]
          mass_to_comp <- recap_tab[toReturn[[i]][prob_edges][1],"mz"]
          if(all(is.na(mass_to_comp))) next
          
          pmin <- which.min(abs(dmz-mass_to_comp))
          if(length(pmin)!=0){
            to_remove[nrm:(nrm+length(posf)-2)] <- posf[-pmin]
          }
          nrm <- nrm + length(posf)-1
        }
      }

      to_remove <- to_remove[1:(nrm-1)]
      if(nrm>1){
        pf <- pf[-to_remove,]
      }
      if(!is.null(nrow(pf)))
      {
        dag <- add_edges(dag,as.numeric(t(pf)),lab=toReturn[[i]][pf])
      }
      
      ###We increment the counter value.
      if(ecount(dag)>0){
        tlab <- table(edge_attr(dag,"lab"))
        ntlab <- as.numeric(names(tlab))
        recap_tab$count[ntlab] <- recap_tab$count[ntlab]+as.numeric(tlab)
      }
    }
    toReturn[[i]] <- dag
  }
  return(list(spectra=list_spec,dags = toReturn, elems = recap_tab))
}

###Return an extended spectra and possibly add the precurosr if needed.
###If maximum is not NULL return the top "k" or top "k-1" depending of the inclusion
###of the precursor.
#' Return the spectrum considered by mineMS2 when building the graph.
#'
#' @param spec A "Spectrum2" object as defined in MSnBase
#' @param ppm A ppm parameter to match the precursor
#' @param dmz A minimum deviation in dalton to match the precursor.
#' @param removeSup A boolean indicating if peaks with a mass higher than the precursor should be remode
#' @param addPrecursor A boolen indicating if a preucrsor should be added if needed
#' @param maxNum The maximum number of peak to keep. The value correspond ot the top-k intensities.
#' @param relInt Should the relative intensity be added to the dataset.
#'
#' @return A data.frame containg the mz, intensity and precursor columns.
#' @export
#'
#' @examples
#' #Loading data
#' data(m2l)
#' 
#' getExtendedSpec(m2l["S10"],maxNum = 15)
getExtendedSpec <- function(spec,ppm=5,dmz=0.02,removeSup=TRUE,addPrecursor=TRUE, maxNum = NULL, relInt = TRUE){
  
  
  ###Reconstruct a spectra form an msnbase object.
  df_spec <- data.frame(mz=mz(spec),intensity=intensity(spec))
  
  if(nrow(df_spec)==0){
    return(data.frame(mz=numeric(0),intensity=numeric(0),precursor=logical(0)))
  }
  
  mzprec <- precursorMz(spec)
  tol <- max(mzprec*ppm*(10e-6),dmz)
  
  ###We remove the values with a masses higher than the precursor.
  if(removeSup){
    psup <- which(df_spec$mz<(mzprec+tol))
    if(length(psup)>0) df_spec <- df_spec[psup,,drop=FALSE]
  }
  
  ###intensity is the minimal intensity of the spectra if the peak does not exists.
  intprec <- min(df_spec$int)
  
  ###Checking if there is a potential precursor
  vec_prec <- NULL
  pprec <- which(abs(df_spec$mz-mzprec)<tol)
  if(length(pprec)==0){
    vec_prec <-rep(FALSE,nrow(df_spec))
    if(addPrecursor){
      df_spec[nrow(df_spec)+1,] <- c(mzprec,intprec)
      vec_prec[nrow(df_spec)] <- TRUE
    }
  }else{
    vec_prec <- rep(FALSE,nrow(df_spec))
    vec_prec[pprec[which.min(df_spec$mz[pprec]-mzprec)]]<-TRUE
  }
  df_spec$precursor <- vec_prec
  pos_prec <- which(vec_prec)
  
  ###Precursor is included in every case.
  if(!is.null(maxNum)){
    if(maxNum < nrow(df_spec)){
      pm <- order(df_spec$int,decreasing=TRUE)
      tpos <- pm[1:maxNum]
      if(addPrecursor & !(pos_prec %in% tpos)){
        tpos <- c(tpos[1:(maxNum-1)],pos_prec)
      }
      df_spec <- df_spec[tpos,]
    }
  }
  
  if(relInt){
    onames <- colnames(df_spec)
    onames <- c(onames,"rel_int")
    rel_int <- (as.numeric(df_spec$int)/sum(as.numeric(df_spec$int)))*100
    # if(any(is.na(rel_int))) browser()
    df_spec$rel_int <- rel_int
  }
  df_spec <- df_spec[order(df_spec$mz,decreasing=TRUE),]
  df_spec
}

###Discretize set of values furnished in the vmz vector
###with a bandwidth given as the minimum between ppm and dmz.
###frac give the minimum number of samples.
discretizeSequenceByClasses <- function(ppm,
                                        vmz, vsample, ndiff, frac = 0,
                                        inter = 0.2, dmz = 0.01,
                                        nPoints = 512, extended=FALSE) {
  porder <- order(vmz)
  vall <- vmz[porder]
  mzrange <- range(vall)
  mzrange <-
    c(max(floor(mzrange[1] / 10) * 10, 1), ceiling(mzrange[2] / 10) * 10)
  massInter <- seq(mzrange[1], mzrange[2], inter / 2)
  masspos <- FindEqualGreaterM(vall, massInter)
  num_group <- 0
  idxgroup <- vector(mode = "list", 512)
  listgroup <- matrix(0, nrow = 512, ncol = 4)
  pos = 0
  previousMz = NULL
  num_group = 0
  for (i in seq(length = length(massInter) - 2)) {
    mem = i
    start <- masspos[i]
    end <- masspos[i + 2] - 1
    if (end - start < 0)
      next
    subpeakl = vall
    pos = pos + 1
    bw = max((ppm * massInter[i] * 1e-6), dmz)/3
    den <-
      density(
        subpeakl,
        bw,
        from = massInter[i] - 3 * bw,
        to = massInter[i + 2] +
          3 * bw,
        n = nPoints
      )
    
    
    bmin1 <- massInter[i] - 3 * bw
    bmax1 <- massInter[i+2] + 3 * bw
    
    ###Putting the close to 0 to 0
    maxy = max(den$y)
    den$y[which(den$y <= bw * maxy)] = 0
    plim = c(-1, 0, 0)
    if (den$y[2] > den$y[1])
      plim[3] <- 1
    
    
    
    oldMz = previousMz
    repeat {
      plim = findLimDensity(den$y, plim[2] + 1, plim[3])
      if (plim[2] <= plim[1])
        break
      selectedPGroup = which(subpeakl >= den$x[plim[1]] &
                               subpeakl <= den$x[plim[2]])
      if (length(selectedPGroup) == 0)
        next
      num_group = num_group + 1
      if (num_group > nrow(listgroup)) {
        listgroup <- rbind(listgroup,
                           matrix(
                             0,
                             nrow = nrow(listgroup),
                             ncol = 4
                           ))
      }
      if (mean(vmz[porder[selectedPGroup]]) %in% oldMz ||
          (abs(median(subpeakl[selectedPGroup]) - massInter[i]) < bw *
           2) ||
          (abs(median(subpeakl[selectedPGroup]) - massInter[i + 2]) <
           bw * 2)) {
        num_group = num_group - 1
        next
      }
      if (length(unique(vsample[porder[selectedPGroup]])) / ndiff < frac) {
        num_group = num_group - 1
        next
      }
      if(extended){
        listgroup[num_group, ] <-
          c(mean(vmz[porder[selectedPGroup]]),
            den$x[plim[1]],
            den$x[plim[2]],bw)
      }else{
        listgroup[num_group, ] <-
          c(mean(vmz[porder[selectedPGroup]]),
            range(vmz[porder[selectedPGroup]]),bw)
      }
      idxgroup[[num_group]] <- porder[selectedPGroup]
      previousMz <- c(previousMz, mean(vmz[porder[selectedPGroup]]))
    }
    
  }
  if (num_group < 1) {
    return(list(dicrete_elem = numeric(0), idx = list()))
  }
  listgroup <- listgroup[1:num_group,]
  idxgroup <- idxgroup[1:num_group]
  ###Reordering by increasing masses.
  og <- order(listgroup[,1],decreasing=FALSE)
  listgroup <- listgroup[og,]
  idxgroup <- idxgroup[og]
  return(list(discrete_elem = listgroup, idx = idxgroup))
}

##Generate a matrix of differences given a sequence.
generateMatDiff <- function(seqmz){
  rmat <- outer(seqmz, seqmz, FUN = "-")
  rmat[lower.tri(rmat, diag = TRUE)] <- NA
  return(list(
    mz = rmat[upper.tri(rmat)],
    idx = which(upper.tri(rmat), arr.ind = TRUE)
  ))
}


#convertMatrixIgraph <- function(x){
#  x[which(!is.na(x),arr.ind = TRUE)] <- 1
#  x[which(is.na(x),arr.ind = TRUE)] <- 0
#  x
#}


checkFracParam <- function (num, maxVal, OneEqualAll = TRUE)
{
  if (num != round(num)) {
    if (num > 1) {
      stop("Param too high.")
    }
    else {
      return(floor(maxVal * num))
    }
  }
  else {
    if (num > 1 & num <= maxVal)
      return(num)
    if (num > maxVal)
      stop("Param too high.")
    if (num == 1 & OneEqualAll)
      return(maxVal)
    if (num == 1)
      return(1)
    if (num == 0)
      return(0)
  }
}


#sumFormula <- function(f1,f2,all_lab=NULL){
#  if(is.na(f1)|is.na(f2)) return(NA)
#  if(is.null(all_lab)){
#    all_lab <- tabAtoms()$name
#  }
#  vf <- rep(0,length(all_lab))
#  names(vf) <- all_lab
#  vf[all_lab] <- vf[all_lab]+f1
#  vf[all_lab] <- vf[all_lab]+f2
  ###We remover the 0
#  vf <- vf[vf!=0]
#  return(vf)
#}


###Function used for labels fusion if there is an incoherence.
fuseElem <- function(elems,dags,thresh=2,atoms=NULL){
  if(is.null(atoms)) atoms <- tabAtoms()$name
  anyChange <- FALSE
  
  ###We first store a list of all the formula when they exist
  allformula <- sapply(elems$formula,MzDiffFormulaFromSingleString,ref=atoms,sep="|")
  # omzmin <- elems$mzmin
  # omzmax <- elems$mzmax
  
  ###Data structure for triple value, a sparse matrix.
  stormat <- Matrix(0,nrow=nrow(elems),ncol=nrow(elems),sparse = TRUE)
  storval <- vector(mode="list",length=1000)
  nval <- 1
  
  ####We get a list of all the possible errors in the labels.
  for(igl in 1:length(dags)){
    g <- dags[[igl]]
    na_f <- TRUE
    if (!is.null(g))
    {
      for(e in is.na(g))
      {
        if(!e) na_f <- FALSE
      }
    }
    
    if(is.null(g) || (length(g) == 1 && na_f)){
      next
    }
    adj_list <- adjacent_vertices(g,V(g),"out")
    lab_sec <- E(g)$lab
    for(n in V(g)){
      adjv <- adj_list[[n]]
      if(length(adjv)>0){
        for(na in adjv){
          succ <- intersect(adj_list[[n]],adj_list[[na]])
          if(length(succ)!=0){
            for(s in succ){
              all_ids <- get_edge_ids(g,c(n,na,na,s,n,s))
              la <- lab_sec[all_ids[1]]
              lb <- lab_sec[all_ids[2]]
              lc <- lab_sec[all_ids[3]]
              
              ma <- min(la,lb)
              mb <- max(la,lb)
              
              ###We check if the value may be added to the dataset or not.
              is_in <- stormat[ma,mb]
              if(is_in==0){
                ###If necessary the storage vector is pushed up
                if(nval>length(storval)){
                  storval <- c(storval,vector(mode="list",length=length(storval)))
                }
                ###The position is added to the dataset.
                stormat[ma,mb] <- nval
                storval[[nval]] <- c(ma,mb,lc)
                nval <- nval+1
              }else{
                if(!lc %in% storval[[is_in]] ){
                  storval[[is_in]] <- c(storval[[is_in]],lc)
                }
              }
            }
          }
        }
      }
    }
  }
  if(nval>1){
    storval <- storval[1:(nval-1)]
  }else{
    storval <- list()
  }
  #print(storval)
  if(length(storval) == 0) return(list(elems=elems,dags=dags,change=anyChange))
  ## Now we have a data structure with all the value, we check which set needs to be fused.
  vmul <- sapply(storval,function(x){length(unique(x))>3})
  pmistake <- which(vmul)
  
  nreduc <- 0
  if(length(pmistake)>0){
    anyChange <- TRUE
    nreduc <- sapply(storval[pmistake],length)-3
  }
  
  ###We fuse the dataset if necessary (if an error has been found)
  to_rm <- numeric(100)
  idrm <- 1
  
  nlab <- 1:nrow(elems)
  
  merged_val <- numeric(100)
  idm <- 1
  
  for(i in pmistake){
    
    sv <- storval[[i]]
    
    a <- sv[1]
    b <- sv[2]
    cs <- sort(sv[3:length(storval[[i]])])
    rcs <- range(cs)
    
    ####We extend the cd label if a value has been jump.
    cs <- rcs[1]:rcs[2]
    
    ###If any of the value is in the merged vals or in the removed vals we skip this step.
    if(any(cs %in% to_rm[1:(idrm-1)])||any(cs %in% merged_val[1:(idm-1)])){
      next
    }
    
    telems <- elems[cs,,drop=FALSE]
    
    ##The dataset are fused.
    mzmin <- min(telems$mzmin)
    mzmax <- max(telems$mzmax)
    mz <- mean(telems$mz)
    fused <- TRUE
    count <- sum(telems$count)
    
    ###We determine correct formula.
    tform <- allformula[cs]
    
    
    
    ####This formula check is not necessary anymore.
    formula <- NULL
    with_form <- which(sapply(tform,length)!=0)
    if(length(with_form)==0){
      formula <- HIGH_MASS_SYMBOL
    }else{
      ###We calculate the new formula.
      tempf <- Reduce(f = addFormula,tform)
      formula <- paste(as.character(tempf),collapse = "|")
    }
    
    ###
    n_add <- length(cs)-1
    if((idrm+n_add-1)>length(to_rm)){
      to_rm <- c(to_rm,numeric(length(to_rm)))
    }
    to_rm[idrm:(idrm+n_add-1)] <- cs[2:length(cs)]
    idrm <- idrm + n_add
    
    if((idm)>length(merged_val)){
      merged_val <- c(merged_val,numeric(length(merged_val)))
    }
    merged_val[idm] <- cs[1]
    
    idm <- idm+1
    elems$mzmin[cs[1]] <- mzmin
    elems$mzmax[cs[1]] <- mzmax
    elems$mz[cs[1]] <- mz
    elems$fused[cs[1]] <- TRUE
    elems$count[cs[1]] <- count
    elems$formula[cs[1]] <- formula
    nlab[cs] <- cs[1]
  }
  
  ###The conserved labels
  ulabs <- sort(unique(nlab))
  
  ###The new labels.
  new_labs <- match(nlab,ulabs)
  
  ###Relabeling of the graph
  for(igl in 1:length(dags)){
    g <- dags[[igl]]
    na_f <- TRUE
    if (!is.null(g))
    {
      for(e in is.na(g))
      {
        if(!e) na_f <- FALSE
      }
    }
    if(is.null(g)||(length(g) == 1 && na_f)){
      next
    }
    if(ecount(g)>0){
      eg_lab <- edge_attr(g,"lab")
      g <- set_edge_attr(g,name = "lab",value = new_labs[nlab[eg_lab]])
    }
    dags[[igl]] <- g
  }
  # checkCorrectness(omzmin,omzmax,elems$mzmin,elems$mzmax,nlab,margin=0.0001)
  
  elems <- elems[ulabs, ]
  elems$lab <- 1:length(ulabs)
  
  tempr <- cleanupElems(elems,dags,thresh)
  
  return(list(elems=tempr$elems,dags=tempr$dags,change=anyChange))
}



cleanupElems <- function(elems,dags,thresh){
  countv <- rep(0,nrow(elems))
  for(igl in 1:length(dags)){
    g <- dags[[igl]]
    
    na_f <- TRUE
    if (!is.null(g))
    {
      for(e in is.na(g))
      {
        if(!e) na_f <- FALSE
      }
    }
    if(is.null(g)||(length(g) == 1 && na_f)){
      next
    }
    eg_lab <- edge_attr(g,"lab")
    if(ecount(g)>0){
      tlab <- table(eg_lab)
      ntlab <- as.numeric(names(tlab))
      countv[ntlab] <- countv[ntlab]+as.numeric(tlab)
    }
  }
  
  p_rm <- elems$lab[countv<thresh]
  p_keep <- elems$lab[countv>=thresh]
  
  nlab <- elems$lab
  nlab[p_keep] <- 1:length(p_keep)
  
  ###Now removing the non frequent edges.
  for(igl in 1:length(dags)){
    g <- dags[[igl]]
    na_f <- TRUE
    if (!is.null(g))
    {
      for(e in is.na(g))
      {
        if(!e) na_f <- FALSE
      }
    }
    
    if(is.null(g)||(length(g) == 1 && na_f)){
      next
    }
    eg_lab <- edge_attr(g,"lab")
    v_rm <- eg_lab %in% p_rm
    if(sum(v_rm)>0){
      g <- delete_edges(g,E(g)[v_rm])
    }
    eg_lab <- edge_attr(g,"lab")
    edge_attr(g,"lab") <- nlab[eg_lab]
    dags[[igl]] <- g
  }
  elems$count <- countv
  elems <- elems[p_keep,]
  return(list(elems=elems,dags=dags))
}

