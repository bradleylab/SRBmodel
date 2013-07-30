plotlines=1
plot(d34SH2S[plotlines,],D3xSH2S[plotlines,],type='l', xlab=xlab.name, ylab=ylab.name, xlim=c(x.axis.min,x.axis.max),ylim=c(y.axis.min,y.axis.max), col='grey30')
plotlines=seq(10,90,9)
for (pcount in plotlines){
lines(d34SH2S[pcount,],D3xSH2S[pcount,], col='grey30')	
}
pcount=100
lines(d34SH2S[pcount,],D3xSH2S[pcount,], col='grey30')

plotlines=c(1,seq(10,100,10))
for (pcount in plotlines){
lines(d34SH2S[,pcount],D3xSH2S[,pcount], col='grey30')	
}

#########
#plot outline only

plotlines=1
plot(d34SH2S[plotlines,],D3xSH2S[plotlines,],type='l', xlab=xlab.name, ylab=ylab.name, xlim=c(x.axis.min,x.axis.max),ylim=c(y.axis.min,y.axis.max), lwd=3, col='grey0', lty=1)
plotlines=100
lines(d34SH2S[pcount,],D3xSH2S[pcount,], lwd=3, col='grey0', lty=1)
lines(d34SH2S[,pcount],D3xSH2S[,pcount], col='grey0', lwd=3, lty=1)	