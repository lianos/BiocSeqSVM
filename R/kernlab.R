SpectrumFeatures <- function(x, length=c(4,5)) {
  stopifnot(inherits(x, "DNAStringSet"))

  if (length(length) == 1) {
    x <- oligonucleotideFrequency(x, length)
  } else {
    x <- do.call(cbind, lapply(length, function(len) {
      oligonucleotideFrequency(x, l)
    }))
  }

  x[, colSums(x) > 0]
}

setMethod("ksvm", c(x="DNAStringSet"),
function(x, y=values(x)$class, kernel="spectrum", length=4, lambda=1.1,
         normalized=TRUE, type=NULL, ...) {
  if (is.null(y)) {
    stop("Label vector y is required")
  }
  kernels <- c('spectrum', 'boundrange', 'constant', 'exponential')
  kernel <- match.arg(kernel, kernels)
  kernel <- stringdot(length=length, lambda=lambda, type=kernel,
                      normalized=normalized)

  if (is.null(type)) {
    warning("SVM type is not specified: trying to guess by labels")
    if (is.integer(y) && length(table(y)) == 2L) {
      type <- 'C-svc'
    }
  }

  ksvm(as.list(as.character(x)), y=y, kernel=kernel, type=type, ...)
})


################################################################################
## Overriding kernlab,predict to "catch" newdata==DNAStringSet
setMethod("predict", signature(object = "ksvm"),
function (object, newdata, type = "response", coupler = "minpair")
{
  ## Modification here by SL to catch DNAStringSet object and turn them
  ## into the kernlab-expected list-of-characters
  if (inherits(kernelf(object), "stringkernel") && !missing(newdata)) {
    if (is(newdata, "DNAStringSet")) {
      newdata <- as.list(as.character(newdata))
    } else if (!is(newdata, "kernelMatrix")) {
      stop("Illegal `newdata` for ksvm,stringkernel")
    }
  }
  ## End modification

  type <- match.arg(type,c("response","probabilities","votes","decision"))
  if (missing(newdata) && type=="response" & !is.null(fitted(object)))
    return(fitted(object))
  else if(missing(newdata))
    stop("Missing data !")

 if(!is(newdata,"list")){
  if (!is.null(kernlab:::terms(object)) & !is(newdata,"kernelMatrix"))
    {
      if(!is.matrix(newdata))
        newdata <- model.matrix(delete.response(kernlab:::terms(object)), as.data.frame(newdata), na.action = n.action(object))
    }
  else
    newdata  <- if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)


    newnrows <- nrow(newdata)
    newncols <- ncol(newdata)
    if(!is(newdata,"kernelMatrix") && !is.null(xmatrix(object))){
      if(is(xmatrix(object),"list") && is(xmatrix(object)[[1]],"matrix")) oldco <- ncol(xmatrix(object)[[1]])
      if(is(xmatrix(object),"matrix")) oldco <- ncol(xmatrix(object))
      if (oldco != newncols) stop ("test vector does not match model !")
    }
  }
  else
   newnrows <- length(newdata)

  p <- 0

  if (is.list(scaling(object)))
    newdata[,scaling(object)$scaled] <-
      scale(newdata[,scaling(object)$scaled, drop = FALSE],
            center = scaling(object)$x.scale$"scaled:center", scale  = scaling(object)$x.scale$"scaled:scale")

  if(type == "response" || type =="decision" || type=="votes")
    {
      if(kernlab:::type(object)=="C-svc"||kernlab:::type(object)=="nu-svc"||kernlab:::type(object)=="C-bsvc")
        {
          predres <- 1:newnrows
          if(type=="decision")
	   votematrix <- matrix(0,kernlab:::nclass(object)*(kernlab:::nclass(object)-1)/2,newnrows)
	  else
	   votematrix <- matrix(0,kernlab:::nclass(object),newnrows)

	  for(i in 1:(kernlab:::nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:kernlab:::nclass(object))
                {
                  p <- p+1

                  if(is(newdata,"kernelMatrix"))
                    ret <- newdata[,which(SVindex(object)%in%alphaindex(object)[[p]]), drop=FALSE] %*% coef(object)[[p]] - b(object)[p]
                  else
                    ret <- kernelMult(kernelf(object),newdata,xmatrix(object)[[p]],coef(object)[[p]]) - b(object)[p]

                  if(type=="decision")
                    votematrix[p,] <- ret
                  else{
                    votematrix[i,ret<0] <- votematrix[i,ret<0] + 1
                    votematrix[j,ret>0] <- votematrix[j,ret>0] + 1
                  }
                }
            }
          if(type == "decision")
            predres <-  t(votematrix)
          else
            predres <- sapply(predres, function(x) which.max(votematrix[,x]))
        }

  if(kernlab:::type(object) == "spoc-svc")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,kernlab:::nclass(object),newnrows)
      for(i in 1:kernlab:::nclass(object)){
        if(is(newdata,"kernelMatrix"))
          votematrix[i,] <- newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% coef(object)[[i]]
        else if (is(newdata,"list"))
          votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],coef(object)[[i]])
        else
          votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],coef(object)[[i]])
      }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(kernlab:::type(object) == "kbb-svc")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,kernlab:::nclass(object),newnrows)
      A <- rowSums(alpha(object))

      for(i in 1:kernlab:::nclass(object))
        {
          for(k in (1:i)[-i])
            if(is(newdata,"kernelMatrix"))
              votematrix[k,] <- votematrix[k,] - (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% alpha(object)[,k][alphaindex(object)[[i]]] + sum(alpha(object)[,k][alphaindex(object)[[i]]]))
            else if (is(newdata,"list"))
              votematrix[k,] <- votematrix[k,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],alpha(object)[,k][alphaindex(object)[[i]]]) + sum(alpha(object)[,k][alphaindex(object)[[i]]]))
            else
              votematrix[k,] <- votematrix[k,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],alpha(object)[,k][alphaindex(object)[[i]]]) + sum(alpha(object)[,k][alphaindex(object)[[i]]]))

          if(is(newdata,"kernelMatrix"))
            votematrix[i,] <- votematrix[i,] + (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% A[alphaindex(object)[[i]]] + sum(A[alphaindex(object)[[i]]]))
          else if (is(newdata,"list"))
            votematrix[i,] <- votematrix[i,] + (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],A[alphaindex(object)[[i]]]) + sum(A[alphaindex(object)[[i]]]))
          else
            votematrix[i,] <- votematrix[i,] + (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],A[alphaindex(object)[[i]]]) + sum(A[alphaindex(object)[[i]]]))

          if(i <= (kernlab:::nclass(object)-1))
            for(kk in i:(kernlab:::nclass(object)-1))
              if(is(newdata,"kernelMatrix"))
                votematrix[kk+1,] <- votematrix[kk+1,] - (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% alpha(object)[,kk][alphaindex(object)[[i]]] + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
              else if (is(newdata,"list"))
                votematrix[kk+1,] <- votematrix[kk+1,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],alpha(object)[,kk][alphaindex(object)[[i]]]) + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
              else
                votematrix[kk+1,] <- votematrix[kk+1,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],alpha(object)[,kk][alphaindex(object)[[i]]]) + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
        }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }
}

  if(type == "probabilities")
    {
      if(is.null(prob.model(object)[[1]]))
        stop("ksvm object contains no probability model. Make sure you set the paramater prob.model in ksvm during training.")

      if(kernlab:::type(object)=="C-svc"||kernlab:::type(object)=="nu-svc"||kernlab:::type(object)=="C-bsvc")
        {
          binprob <- matrix(0, newnrows, kernlab:::nclass(object)*(kernlab:::nclass(object) - 1)/2)
          for(i in 1:(kernlab:::nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:kernlab:::nclass(object))
                {
                  p <- p+1
                  if(is(newdata,"kernelMatrix"))
                    binprob[,p] <- 1 - .SigmoidPredict(as.vector(newdata[,which(SVindex(object)%in%alphaindex(object)[[p]]), drop=FALSE] %*% coef(object)[[p]] - b(object)[p]), prob.model(object)[[p]]$A, prob.model(object)[[p]]$B)
                  else
                    binprob[,p] <- 1 - .SigmoidPredict(as.vector(kernelMult(kernelf(object),newdata,xmatrix(object)[[p]],coef(object)[[p]]) - b(object)[p]), prob.model(object)[[p]]$A, prob.model(object)[[p]]$B)
                }
            }
          multiprob <- couple(binprob, coupler = coupler)
        }
      else
        stop("probability estimates only supported for C-svc, C-bsvc and nu-svc")
    }

  if(kernlab:::type(object) == "one-svc")
    {
      if(is(newdata,"kernelMatrix"))
        ret <- newdata %*% coef(object) - b(object)
      else
        ret <- kernelMult(kernelf(object),newdata,xmatrix(object),coef(object)) - b(object)
      ##one-class-classification: return TRUE/FALSE (probabilities ?)
      if(type=="decision")
      	return(ret)
	else
	{
	ret[ret>0]<-1
      	return(ret == 1)
      }
    }
  else {
    if(kernlab:::type(object)=="eps-svr"||kernlab:::type(object)=="nu-svr"||kernlab:::type(object)=="eps-bsvr")
      {
        if(is(newdata,"kernelMatrix"))
          predres <- newdata %*% coef(object) - b(object)
       else
         predres <- kernelMult(kernelf(object),newdata,xmatrix(object),coef(object)) - b(object)
      }
    else {
      ##classification & votes : return votematrix
      if(type == "votes")
        return(votematrix)

      ##classification & probabilities : return probability matrix
      if(type == "probabilities")
        {
          colnames(multiprob) <- lev(object)
          return(multiprob)
        }

      if(is.numeric(lev(object)) && type == "response")
         return(lev(object)[predres])

      if (is.character(lev(object)) && type!="decision")
        {
          ##classification & type response: return factors
          if(type == "response")
            return(factor (lev(object)[predres], levels = lev(object)))
        }
    }
  }

  if (!is.null(scaling(object)$y.scale) & !is(newdata,"kernelMatrix") & !is(newdata,"list"))
    ## return raw values, possibly scaled back
    return(predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center")
  else
    ##else: return raw values
    return(predres)
})
