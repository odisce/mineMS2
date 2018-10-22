

###Convert components furnished to spectral ids.
convertComponent <- function(comp,ids){
	if(is.numeric(comp)){
		trad <- ids[comp]
		if(any(is.na(trad))) stop("Invalid component elements: ",comp[is.na(trad)])
	}else{
		trad <- match(comp,ids)
		if(any(is.na(trad))) stop("Invalid component elements: ",comp[is.na(trad)])
		return(paste("S",trad,sep=""))
	}
}


#' Find pattern explainig components
#'
#' Find pattern maximizing the socre givne in metrics with each components. Threshold may be used to indicate that not match
#' has been found. The majority
#'
#' @param m2l An ms2Lib object
#' @param components A list containing the components to be explained.
#' @param metric The metric used to match a component and the pattern occurence.
#' @param threshold The threshold used to discriminate the dataset
#'
#' @return A list containing one data.frame by component including the following informations:
#'#' \description{
#' \item[id] The id of the matched patterns.
#' \item[f1] The F1-score between the pattern and the component.
#' \item[accuracy] The accuracy between the pattern and the component.
#' \item[recall] The recall between the pattern and the component.
#' \item[miss] The miss rate between the pattern and the component.
#' }
#' @export
#'
#' @examples
#' print("Examples to be put here")
findPatternsFromComponents <- function(m2l,components,metric=c("f1"),threshold=NA_real_){


	###We first convert the ids of the compoentns to creect format.
	components <- sapply(components,convertComponent,ids=mm2Ids(m2l),simplify = FALSE)

	###We determine the best format for each found pattern.

	bpat <- sapply(component,function(x,m2l,crit){
		find.patterns.class(m2l,x,criterion = crit,full = FALSE)
	},m2l=m2l,crit=metric,simplify = FALSE)

	###Postprocessing
	return(bpat)
}