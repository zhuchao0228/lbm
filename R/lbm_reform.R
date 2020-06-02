lbm.reform<-function(x,y, t = NULL, num_bvectors = NULL,
                    indicator = NULL,start,vce,control,...)
{
  # sub_x<-subset(x,select=-`(Intercept)`)
  sub_x<-subset.matrix(x,select=-1)
  # print(sub_x)
  rownames.sub_x<-row.names(sub_x)
  # t <- as.data.frame(bvectors)
  rownames.bv<-row.names(t)
  colnames.bv<-colnames(t)
  t.repa<-as.data.frame(matrix(NA,ncol=ncol(t),
                        nrow=nrow(t),byrow=T,
                        dimnames=dimnames(t)))
  for(r in 1:num_bvectors)
  {
    if(r == 1) {
      sub_x<-matrix(unlist(apply(sub_x,1,function(p) p-t[r,])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t.repa[r,]<-temp<-t[r, ]
      t<-matrix(unlist(apply(t,1,function(p) p-temp)),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
    } else {
      q <- r-1
      if(t[r,q]!=0) {
        temp <- t[r, ]/t[r, q]
      } else {
        if(r!=num_bvectors) {
          nonzero<-apply(t[r:num_bvectors,1:ncol(sub_x)],
                        2, function(a) which(a!=0))
          max.nonzero.col<-which.max(unlist(lapply(nonzero,length)))
          nonzero.row<-which(round(t[,max.nonzero.col],digits=15)!=0)
          name.temp<-rownames.bv[r]
          rownames.bv[r]<-rownames.bv[nonzero.row[1]]
          rownames.bv[nonzero.row[1]]<-name.temp
        } else {
          nonzero<-which(t[r:num_bvectors,]!=0)
          max.nonzero.col<-nonzero[1]
        }
        # The column which has the most none zero members involved.
        name.temp<-colnames.bv[q]
        colnames.bv[q]<-colnames.bv[max.nonzero.col]
        colnames.bv[max.nonzero.col]<-name.temp
        sub_x<-sub_x[,colnames.bv]
        t<-t[rownames.bv,colnames.bv]
        t.repa<-t.repa[,colnames.bv]
        temp <- t[r, ]/t[r, q]
      }
      t.repa[r,]<-temp
      sub_x<-matrix(unlist(apply(sub_x,1,function(p) p-temp*p[q])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t<-matrix(unlist(apply(t,1,function(p) p-temp*p[q])),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
      # print(sub_x)
    }
  }
  temp_x <- cbind(indicator=indicator, sub_x)
  temp_x <- temp_x[-which(temp_x[,"indicator"]==1), ]
  y<-y[-which(indicator==1)]
  sub_x<-subset(temp_x,select=-c(1:num_bvectors))
  results_temp<-tryCatch({
      res<-lbm.NR(x=sub_x,y=y, vce=vce,control=control,
                  start=start[colnames(sub_x)],intercept=FALSE,...)
    },
      error=function(err){
        return(ind=0)
    })
  if(is.list(results_temp)) {
    vcov<-results_temp$vcov
    beta <- results_temp$coefficients
    res <- list(vcov = vcov, beta = beta, t.repa=t.repa,
                deviance=results_temp$deviance,fitted.values=results_temp$fitted.values)
    # print(1)
    return(res)
  } else {
    temp_init<-lbm.InitValue(x=sub_x,y=y,noconstant=T)
    if(is.list(temp_init)) {
      if(temp_init$SAS==0) {
        results_temp<-lbm.NR(x=sub_x,y=y,mu=temp_init$mustart,
                          vce=vce,control=control,intercept=FALSE,...)
      } else {
        results_temp<-lbm.NR(x=sub_x,y=y,start=temp_init$start,
                          vce=vce,control=control,intercept=FALSE,...)
      }
      vcov <- results_temp$vcov
      beta <- results_temp$coefficients
      res <- list(vcov = vcov, beta = beta,t.repa=t.repa,
                  deviance=results_temp$deviance,fitted.values=results_temp$fitted.values)
      # print(2)
      return(res)
    } else {
      return(temp_init)
    }
  }
}
