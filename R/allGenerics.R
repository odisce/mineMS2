###Getter and setter (internal)

setGeneric("mm2Spectra", function(m2l) standardGeneric("mm2Spectra"))

setGeneric("mm2Dags", function(m2l) standardGeneric("mm2Dags"))

setGeneric("mm2SpectraInfos", function(m2l) standardGeneric("mm2SpectraInfos"))

setGeneric("mm2EdgesLabels", function(m2l) standardGeneric("mm2EdgesLabels"))

setGeneric("mm2NodesLabels", function(m2l) standardGeneric("mm2NodesLabels"))


setGeneric("mm2Spectra<-", function(m2l,value) standardGeneric("mm2Spectra<-"))

setGeneric("mm2Dags<-", function(m2l,value) standardGeneric("mm2Dags<-"))

setGeneric("mm2SpectraInfos<-", function(m2l,value) standardGeneric("mm2SpectraInfos<-"))

setGeneric("mm2EdgesLabels<-", function(m2l,value) standardGeneric("mm2EdgesLabels<-"))

setGeneric("mm2NodesLabels<-", function(m2l,value) standardGeneric("mm2NodesLabels<-"))


###Analysis function
setGeneric("discretizeMassLosses", function(m2l,...) standardGeneric("discretizeMassLosses"))
setGeneric("discretizeFragments", function(m2l,...) standardGeneric("discretizeFragments"))