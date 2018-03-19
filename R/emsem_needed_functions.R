#' A product function of each matrix from an array and its transpose
#'
#' This function allows to compute in the EM algorithm the
#' product of each matrix from an array and its transpose.
#' @param L an array.
#' @return an array composed by each computed matrices.
#' @export
#' @examples
#' trans_tab()

trans_tab<- function(L){ 
  if (class(L) != "array" ) cat( "Warning: the argument is not an array.")
  L_new = array(rep(0,(dim(L)[1])*(dim(L)[1])*(dim(L)[3])), dim = c(dim(L)[1], dim(L)[1], dim(L)[3]))
  for(i in 1:dim(L)[3]){
    L_new[,,i] = L[,,i]%*%t(L[,,i])
  }
  return(L_new)
}

#' A product function of a matrix from an array and another matrix
#'
#' This function allows to compute in the EM algorithm the
#' product of each matrix from an array and another matrix.
#' @param L an array.
#' @param DD a matrix.
#' @return an array composed by each computed matrices.
#' @export
#' @examples
#' prod_tab_mat()
prod_tab_mat<- function(DD,L){  
  if (class(L) != "array" ) cat( "Warning: the argument is not an array.")
  if (class(DD) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  L_new = array(rep(0,(dim(DD)[1])*(dim(L)[2])*(dim(L)[3])), dim = c(dim(DD)[1], dim(L)[2], dim(L)[3]))
  for(i in 1:dim(L)[3]){
    L_new[,,i] = DD%*%L[,,i]
  }
  return(L_new)
}

#' A mean function of a product between each element from an array and the transpose of each row from a matrice
#' reconstructed factor from an array and another matrix
#'
#' This function allows to compute in the EM algorithm 
#' the mean on the statistical units of the products between each element of the reconstructed factors and 
#' the corresponding instance in a bloc of observed data.
#' @param Ftilde an array.
#' @param X a matrix.
#' @return a mean value.
#' @export
#' @examples
#' moy_TildeData()
moy_TildeData<- function(Ftilde,X){ 
  if (class(Ftilde) != "array" ) cat( "Warning: the argument is not an array.")
  if (class(X) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  Ftilde_X = array(rep(0),dim=c(dim(Ftilde)[1],dim(X)[2],n))
  for(i in 1:n){
    Ftilde_X[,,i] =  Ftilde[,,i]%*%t(X[i,])
  }
  moy = apply(Ftilde_X, 1:2, mean)
  return(moy)
}

#' A mean function of a product between each row of a matrice and the transpose of another matrice
#'
#' This function allows to compute in the EM algorithm 
#' the mean on the statistical units of the products between each 
#' row of a bloc of observed data and a row of an extra-covariate matrix.
#' @param YX a matrix.
#' @param T a matrix.
#' @return a mean value.
#' @export
#' @examples
#' moy_T_Prim_Data()
moy_T_Prim_Data<- function(YX,T){ 
  if (class(YX) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  if (class(T) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  YX_Tprim = array(rep(0), dim= c(ncol(YX) , ncol(T) , nrow(YX) ))
  for(i in 1:nrow(YX)){
    YX_Tprim[,,i] =  t(t(YX[i,]))%*%T[i,]#t(t()) allows to obtain a matrix object and not a vector.
  }
  moy = apply(YX_Tprim, 1:2, mean)
  return(moy)
}

#' A mean function of a product between each row of a matrice and its transpose
#'
#' This function allows to compute in the EM algorithm 
#' the mean on the statistical units of the products between each 
#' row of an extra-covariate matrix and a row of the same extra-covariate matrix.
#' @param T a matrix.
#' @param Tbis a matrix.
#' @return a mean value.
#' @export
#' @examples
#' moy_T_T_Prim()
moy_T_T_Prim<- function(T,Tbis){ 
  if (class(T) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  if (class(Tbis) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  T_T_prim = array(rep(0), dim= c(ncol(T) , ncol(T) , nrow(T) ))
  for(i in 1:nrow(T)){
    T_T_prim[,,i] =  t(t(T[i,]))%*%Tbis[i,]
  }
  moy = apply(T_T_prim, 1:2, mean)
  return(moy)
}

#' A mean function of a product between each element from an array and the transpose from another array
#'
#' This function allows to compute in the EM algorithm 
#' the mean on the statistical units of the products between each 
#' element of an explanory reconstructed factor and each element from 
#' a dependent reconstructed factor.
#' @param Ftilde an array corresponding to the explanatory factors.
#' @param Gtilde an array corresponding to the dependent factor.
#' @return a mean value.
#' @export
#' @examples
#' moy_fg_tilde()
moy_fg_tilde <- function(Ftilde, Gtilde){
  if (class(Ftilde) != "array" ) cat( "Warning: the argument is not an array.")
  if (class(Gtilde) != "array" ) cat( "Warning: the argument is not an array.")
  FGtilde = array(rep(0), dim = c(dim(Ftilde)[1], dim(Gtilde)[1], n))
  for(i in 1:n){
    FGtilde[,,i] = Ftilde[,,i]%*%t(Gtilde[,,i])
  }
  moy = apply(FGtilde, 1:2, mean)
  return(moy)
}

#' A product function of each element from an array and each row of a matrix
#'
#' This function allows to compute in the EM algorithm 
#' the products between each matrix from an array and 
#' each row of a matrix.
#' @param mat a matrix.
#' @param L an array.
#' @return an array.
#' @export
#' @examples
#' prod_tab_rowmat()
prod_tab_rowmat<- function(mat,L){ 
  if (class(L) != "array" ) cat( "Warning: the argument is not an array.")
  if (dim(L)[3] != nrow(mat) ) cat( "Warning: the number of matrices from L is different from the number of rows.")
  if (class(mat) != "matrix" ) cat( "Warning: the argument is not a matrix.")
  L_new = matrix(rep(0,(ncol(mat)*nrow(mat))), ncol = ncol(mat), nrow = nrow(mat))
  for(i in 1:nrow(mat)){
    L_new[i,] = mat[i,]*L[,,i]
  }
  return(L_new)
}

#' A scale and center function adapted to vectors
#'
#' This function allows to scale and center a vector.
#' @param v_test a vector.
#' @return a vector.
#' @export
#' @examples
#' cr()
cr<- function(v_test){ 
  if (class(v_test) != "numeric" ) cat( "Warning: the argument is not a vector and you need a numeric class.")
  esp<-mean(v_test)
  ectyp<-sqrt(var(v_test)*(length(v_test)-1)/length(v_test))
  v_test_centre_redui <- (v_test-rep(esp,length(v_test)))/rep(ectyp, length(v_test))
  return(v_test_centre_redui)
}
