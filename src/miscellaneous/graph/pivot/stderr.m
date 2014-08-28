function err=stderr(x)
err=nanstd(x)/sqrt(size(x,1));