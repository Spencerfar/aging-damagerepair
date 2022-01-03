midpoints <- function(x, dp=2){
    lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
    upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
    return(round(lower+(upper-lower)/2, dp))
}


se <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))


invsoftplus <- function(x){
    return(log(exp(x)-1))
}

