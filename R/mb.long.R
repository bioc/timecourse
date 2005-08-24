## Main function for 1-sample and 2-sample longitudinal microarray time course data analysis
## Yu Chuan Tai
## University of California, Berkeley
## Created: August 2004
## Last Modified: August 14, 2005


mb.long <- function(object, method=c("1D", "paired", "2D"), type=c("none", "robust"), times, 
reps, prior.df=NULL, prior.COV=NULL, prior.eta=NULL,
condition.grp=NULL, rep.grp=NULL, time.grp=NULL, one.sample=FALSE, ref=NULL, p=0.02, out.t=FALSE, 
tuning=1.345, HotellingT2.only=TRUE)
{

   
     method <- match.arg(method, c("1D","paired","2D"))
     type <-  match.arg(type, c("none", "robust"))

     if(method=="1D" & type=="none")  out <- mb.1D(object, times, reps, prior.df, prior.COV, 
     prior.eta, time.grp, rep.grp, vec=out.t, prop=p, T2.only=HotellingT2.only)
     if(method=="1D" & type=="robust")  out <- mb.1D(object, times, reps, prior.df, 
     prior.COV, prior.eta, time.grp, rep.grp, r=TRUE, vec=out.t, d=tuning, prop=p, 
     T2.only=HotellingT2.only)
     if(method=="paired" & type=="none")  out <- mb.paired(object, times, reps, prior.df, 
     prior.COV, prior.eta,condition.grp, time.grp, rep.grp, one.sample, ref, vec=out.t, 
     prop=p,T2.only=HotellingT2.only)
     if(method=="paired" & type=="robust")  out <- mb.paired(object, times, reps, prior.df, 
     prior.COV, prior.eta, condition.grp, time.grp, rep.grp, one.sample, ref, r=TRUE, vec=out.t, 
     d=tuning, prop=p, T2.only=HotellingT2.only)
     if(method=="2D" & type=="none") out <- mb.2D(object, times, reps, prior.df, 
     prior.COV, prior.eta,condition.grp, time.grp, rep.grp, vec=out.t, prop=p, T2.only=HotellingT2.only)
     if(method=="2D" & type=="robust") out <- mb.2D(object, times, reps, prior.df, 
     prior.COV, prior.eta, condition.grp, time.grp, rep.grp, r=TRUE, vec=out.t, d=tuning, prop=p, 
     T2.only=HotellingT2.only)
     
     out
}


