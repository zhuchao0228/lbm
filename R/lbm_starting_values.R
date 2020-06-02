
lbm.InitValue<-function(x,y,noconstant=FALSE)
{
  for(u in 1:3)
  {
    if(u==1) {
      ug<-(y+0.5)/2
      z<-log(ug)+(y-ug)/ug
    } else if(u==2) {
      ug<-(y+mean(y))/2
      z<-log(ug)+(y-ug)/ug
    } else {
      ug<-abs(y-0.1)
      z<-log(ug)
    }
    w<-ug/(1-ug)
    init_mu_old<-stats::lm.wfit(x=x,y=z,w=w)
    mu_old<-exp(init_mu_old$fitted.values)
    if(all(mu_old<1 & mu_old>0)) {
      conv<-1
      if(noconstant==FALSE) start<-init_mu_old$coefficients else {
        if(u==3) coe<-init_mu_old$coefficients else mustart<-ug
      }
    } else {
      if(noconstant==FALSE)
      {
        conv<-1
        init_mu_old$coefficients[1]<-init_mu_old$coefficients[1]-log(max(mu_old))+log(0.9)
        start<-init_mu_old$coefficients
        # print(exp(x%*%start))
      } else {
        conv<-0
        #mu_old<-exp(init_mu_old$fitted.values-max(init_mu_old$fitted.values))
        mu_old[which(mu_old>=1)]<-0.9999
        for (p in 1:25)
        {
          if(p==1) {
            ug<-mu_old
          } else {
            ug<-mu_new
          }
          z<-log(ug)+(y-ug)/ug
          w=ug/(1-ug)
          init_mu_new<-stats::lm.wfit(x=x,y=z,w=w)
          mu_new<-exp(init_mu_new$fitted.values)

          if(any(mu_new>1 | mu_new<0)) {
            #mu_new<-exp(init_mu_new$fitted.values-max(init_mu_new$fitted.values))
            mu_new[which(mu_new>=1)]<-0.9999
          } else {
            mustart<-ug
            conv<-1
            break
          }
          if(p==25 & conv!=1) conv<-0
        }
      }
    }
    if(conv==1) break
  }
  if(conv==1)
  {
    if(noconstant==FALSE) return(list(start=start)) else {
      if(u!=3) return(list(mustart=mustart,SAS=0)) else return(list(start=coe,SAS=1))
    }
  } else return(conv=0)
}
