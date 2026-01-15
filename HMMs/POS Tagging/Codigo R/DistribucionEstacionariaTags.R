library(matrixcalc)
library(MASS)

Gamma <- fractions(matrix.power(Gamma, 100))

plot(0,0, ylim=c(0,0.7), xlim=c(0,17))
for (i in 1:17) 
{
  lines(1:17, Gamma[i,], col=i)
}