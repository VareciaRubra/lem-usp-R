library(mnormt)
library(ellipse)


g1<-rmnorm(100, mean=c(7,5), varcov=matrix(c(1,0.5,0.5,1),2,2))
g2<-rmnorm(100, mean=c(3,7), varcov=matrix(c(1,0.5,0.5,1),2,2))


layout(matrix(1:3,1,3))
plot(NA,type="n", xlim=c(0,10), ylim=c(2,10), xlab="Trait 1", ylab="Trait 2")
points(g1, pch=2)
points(g2, pch=3)
polygon(ellipse(var(g1),centre= c(7,5)))
polygon(ellipse(var(g2),centre= c(3,7)))

plot(NA,type="n", xlim=c(-5,5), ylim=c(-5,5), xlab="Trait 1", ylab="Trait 2")
abline(h=0, lty=2)
abline(v=0, lty=2)

polygon(ellipse(var(g1)))
polygon(ellipse(var(g2)))

plot(NA,type="n", xlim=c(-5,5), ylim=c(-5,5), xlab="Trait 1", ylab="Trait 2")
abline(h=0, lty=2)
abline(v=0, lty=2)

lmt<-lm(rbind(g1,g2)~rep(c("A","B"), each=100))

varRes<-(t(lmt$residuals)%*%lmt$residuals)/lmt$df.residual
polygon(ellipse(varRes), lty=3)

dev.copy2pdf(file="fig1.pdf")
dev.print(png, file="fig1.pdf")

###########################


g1<-rmnorm(100, mean=c(9,7), varcov=matrix(c(1,-0.5,-0.5,1),2,2))
g2<-rmnorm(100, mean=c(7,6), varcov=matrix(c(1,0.5,0.5,1),2,2))


layout(matrix(1:3,1,3))
plot(NA,type="n", xlim=c(4,12), ylim=c(2,10), xlab="Trait 1", ylab="Trait 2")
points(g1, pch=2)
points(g2, pch=3)
polygon(ellipse(var(g1),centre= c(9,7)))
polygon(ellipse(var(g2),centre= c(7,6)))

plot(NA,type="n", xlim=c(-5,5), ylim=c(-5,5), xlab="Trait 1", ylab="Trait 2")
abline(h=0, lty=2)
abline(v=0, lty=2)

polygon(ellipse(var(g1)))
polygon(ellipse(var(g2)))

plot(NA,type="n", xlim=c(-5,5), ylim=c(-5,5), xlab="Trait 1", ylab="Trait 2")
abline(h=0, lty=2)
abline(v=0, lty=2)

lmt<-lm(rbind(g1,g2)~rep(c("A","B"), each=100))

varRes<-(t(lmt$residuals)%*%lmt$residuals)/lmt$df.residual
polygon(ellipse(varRes), lty=3)


dev.copy2pdf(file="fig3.pdf")
