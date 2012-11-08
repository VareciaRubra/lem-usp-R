r2.ic<- function(group1, group2, resample.size= 1000)  
{
    group1 <- as.matrix(group1)
    cor.obs1 <- cor(group1, y=NULL, method="pearson")
    r2.obs.1 <- mean((cor.obs1[lower.tri(cor.obs1)])^2)
    r2.est.1 <- numeric(resample.size)
    nr1 <- dim(group1)[1]

    group2<- as.matrix(group2)
    cor.obs2 <- cor(group2, y=NULL, method="pearson")
    r2.obs.2 <- mean((cor.obs2[lower.tri(cor.obs2)])^2)
    r2.est.2 <- numeric(resample.size)
    nr2 <- dim(group2)[1]

    for (i in 1:resample.size)
    {
        group.1<- group1[sample.int(nr1, replace=TRUE),]
        cor.est.1<- cor(group.1, y=NULL, method= "pearson")
        r2.est.1 [i]<- mean(cor.est.1[lower.tri(cor.est.1)]^2) 

        group.2<- group2[sample.int(nr2, replace=TRUE),] 
        cor.est.2<- cor(group.2, y=NULL, method= "pearson")
        r2.est.2 [i]<- mean(cor.est.2[lower.tri(cor.est.2)]^2)
    }

    # Plotting
    par (mar=c(5,4,4,3.5))
    limites<- range(r2.est.1, r2.est.2)
    plot(density(r2.est.1), main="curva de densidade estimada", xlab= "valores de r2", ylab= "densidade", 
         xlim=limites,ylim=c(0,70), col="blue")
    abline(v=r2.obs.1, col="blue")
    quantil.1 <- quantile(r2.est.1, probs= c(0.025,0.975))
    abline(v=quantil.1[1], lty=2, col="blue")
    abline(v=quantil.1[2], lty=2, col="blue")
    par(new=TRUE)
    plot(density(r2.est.2), main="", xlab= "", ylab= "", xlim=limites, ylim=c(0,70), col="red")
    abline(v=r2.obs.2, col="red")
    quantil.2<- quantile(r2.est.2, probs= c(0.025,0.975))
    abline(v=quantil.2[1], lty=2, col="red")
    abline(v=quantil.2[2], lty=2, col="red")
    mtext(text= c("group 1", "\n group 2"),col= c("blue", "red"), side=1, line=1.5, adj=0, padj=1)

    #Results

    vetor.resultados<- list(r2.obs= c(r2.obs.1,r2.obs.2), dist1= c(r2.est.1), 
                            dist2= c(r2.est.2), IC=c(group1= quantil.1, group2= quantil.2))
    return(vetor.resultados)
}
