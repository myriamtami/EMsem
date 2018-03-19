#' An simultaneous EM estimation of a structural equation model and its latent factors.
#'
#' This function allows you to fit a structural equation model with one bloc of dependant observed variables
#' and two blocs of explanatory observed variables. 
#' Underlying to each bloc there is a latent factor which is simultaneously reconstructed.
#' @param Y numeric output matrix (or vector) which is the dependant bloc of observed variables.
#' @param X1 numeric input matrix (or vector) which is one explanatory bloc of observed variables.
#' @param X2 numeric input matrix (or vector) which is one explanatory bloc of observed variables.
#' @param T numeric input matrix (or vector) for extra-covariates linked to Y.
#' @param T1 numeric input matrix (or vector) for extra-covariates linked to X1.
#' @param T2 numeric input matrix (or vector) for extra-covariates linked to X2.
#' @param epsilon numeric accuracy for the convergence criterion.
#' @param nb_it numeric for give a maximal number of iteration.
#' @return list giving the reconstructed latent factor values,
#'  the estimated model’s coefficients, the value of the stopping criterion,
#'  the number of iterations effected.
#' @keywords latent factors
#' @author Myriam Tami, \email{myriam.tami@@univ-grenoble-alpes.fr}
#' @references \url{https://myriamtami.github.io/data/publications/Tami_thesis.pdf}
#' @export
#' @examples
#' emsem_function()

#C'est quoi export?
#je met des mots clés?
#Comment sourcer les fonctions de l'autre .R devant ma fonction?
#Est ce que je met un exemple de datas simulées ?

emsem_function <- function(Y,X1,X2,T,T1,T2,epsilon=10^{-3},nb_it=100){
  
  
  #### Function arguments
  # X1      : (matrix) block of data concerning the first factor 
  # X2      : (matrix) block of data concerning the second factor 
  # Y       : (matrix/vecteur) response variable
  # epsilon : (real) criterion of convergence 
  # nb_it   : (integer) number of iteration
  #nb_it = 100  # limite du nombre d'itirations
  
  ##################################
  #----ENTRE DE PARAMETRES----
  ##################################
  qX1 = ncol(X1)
  qX2 = ncol(X2)
  qY = ncol(Y)
  qT= ncol(T) #New
  qT1= ncol(T1) #New
  qT2= ncol(T2) #New
  kF1 = 1 # CHOIX SELON NOTRE REDACTION : SI != 1 LES RESULTATS NE SONT PAS GARANTIS!
  kF2 = 1
  kG = 1 # CHOIX SELON NOTRE REDACTION : SI != 1 LES RESULTATS NE SONT PAS GARANTIS!
  n = nrow(Y) 
  nb_fact = kF1+kF2+kG
  
  ######################################################
  #----PARAMETRES A ESTIMER : CREATION OBJETS----
  ######################################################
  
  A1 = 0 * matrix(1:(kF1*qX1), nrow = kF1, ncol = qX1)
  A2 = 0 * matrix(1:(kF2*qX2), nrow = kF2, ncol = qX2)
  B = 0 * matrix(1:(kG*qY), nrow = kG, ncol = qY)
  C1 = 0
  C2 = 0
  #New
  D = 0 * matrix(1:(qT*qY), ncol = qY, nrow = qT)
  D1 = 0 * matrix(1:(qT1*qY), ncol = qY, nrow = qT1)
  D2 = 0 * matrix(1:(qT2*qY), ncol = qY, nrow = qT2)
  Psi_X1 = 0 * matrix(1:(qX1*qX1), nrow = qX1, ncol= qX1)
  Psi_X2 = 0 * matrix(1:(qX2*qX2), nrow = qX2, ncol= qX2)
  Psi_Y = 0 * matrix(1:(qY*qY), nrow = qY, ncol = qY)
  sigma2_X1_chap = 0
  sigma2_X2_chap = 0
  sigma2_Y_chap = 0
  
  #Objets necessaires a l'estimation des sigma2 chapeau 
  #Pour X1
  dud1_X1 <- matrix(rep(0, n*qX1), nrow = n, ncol = qX1)
  res1_X1 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res2_X1 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res3_X1 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res_X1_chap <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  #sigma2_X1_chap : objet de taille 1*1 => non construit
  #Pour X2
  dud1_X2 <- matrix(rep(0, n*qX2), nrow = n, ncol = qX2)
  res1_X2 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res2_X2 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res3_X2 <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res_X2_chap <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  #sigma2_X_chap : objet de taille 1*1 => non construit
  #Pour Y
  dud1_Y <- matrix(rep(0, n*qY), nrow = n, ncol = qY)
  res1_Y <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res2_Y <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res3_Y <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  res_Y_chap <- matrix(rep(0, n*1), nrow = n, ncol = 1)
  #sigma2_Y_chap : objet de taille 1*1 => non construit
  
  #Objets necessaires pour le stock des valeurs des parametres a l'iteration precedente : 
  #Objectif : calculer la difference entre chaque it nommee diff pour le critere d'arret
  A1_it = 0 * matrix(1:(kF1*qX1), nrow = kF1, ncol = qX1)
  A2_it = 0 * matrix(1:(kF2*qX2), nrow = kF2, ncol = qX2)
  B_it = 0 * matrix(1:(kG*qY), nrow = kG, ncol = qY)
  C1_it = 0
  C2_it = 0
  #New
  D_it = 0 * matrix(1:(qT*qY), ncol = qY, nrow = qT)
  D1_it = 0 * matrix(1:(qT1*qY), ncol = qY, nrow = qT1)
  D2_it = 0 * matrix(1:(qT2*qY), ncol = qY, nrow = qT2)
  sigma2_X1_chap_it =  0
  sigma2_X2_chap_it =  0
  sigma2_Y_chap_it = 0
  
  
  #####################################################################################################
  #
  #-----INITIALISATION PAR ACP ET REGRESSION----- 
  #
  #####################################################################################################
  #Pour eq de X1
  #Regression pour estimation de D1_sim
  regX1_Fact <-lm(X1 ~ T1 - 1)
  ini_D1 = coef(regX1_Fact)[1:qT1,] #New
  
  
  #Regression pour estimation de A1_sim et sigma2X1
  cl1X1_Fact <- PCA(X1 - T1%*%ini_D1, scale.unit=FALSE)$ind$coord[,1]#donne la premiere composante d'apres FactomineR : correspond au facteur F
  cl1X1_Fact <- cr(cl1X1_Fact)
  regX1_Fact <-lm(X1 - T1%*%ini_D1 ~ cl1X1_Fact - 1 )
  ini_A1_Fact = matrix(coef(regX1_Fact), nrow=1,ncol = qX1, byrow=TRUE)
  sig_lm_qX1_Fact <- rep(NA,qX1)
  for(i in 1:qX1){
    sig_lm_qX1_Fact[i] <- summary(regX1_Fact)[[i]]$sigma
  }
  ini_sigma2X1_Fact = mean(sig_lm_qX1_Fact^2) #initialisation du parametre sigma2X1
  
  
  
  #Pour eq de X2
  #Regression pour estimation de D2_sim
  regX2_Fact <-lm(X2 ~ T2 - 1)
  ini_D2 = coef(regX2_Fact)[1:qT2,] #New
  #Regression pour estimation de A2_sim et sigma2X2
  cl1X2_Fact <- PCA(X2 - T2%*%ini_D2, scale.unit=FALSE)$ind$coord[,1]#donne la premiere composante d'apres FactomineR : correspond au facteur F
  cl1X2_Fact <- cr(cl1X2_Fact)
  regX2_Fact <-lm(X2 - T2%*%ini_D2 ~ cl1X2_Fact - 1 )
  ini_A2_Fact = matrix(coef(regX2_Fact), nrow=1,ncol = qX2, byrow=TRUE)
  sig_lm_qX2_Fact <- rep(NA,qX2)
  for(i in 1:qX2){
    sig_lm_qX2_Fact[i] <- summary(regX2_Fact)[[i]]$sigma
  }
  ini_sigma2X2_Fact = mean(sig_lm_qX2_Fact^2) #initialisation du parametre sigma2X1
  
  #Regression pour estimation de D_sim
  regY_Fact <-lm(Y ~ T - 1)
  ini_D = coef(regY_Fact)[1:qT,] #New 
  #Regression pour estimation de B_sim et sigma2Y
  cl1Y_Fact <- PCA(Y - T%*%ini_D, scale.unit=FALSE)$ind$coord[,1]#donne la premiere composante d'apres FactomineR : correspond au facteur G
  
  cl1Y_Fact <- cl1Y_Fact/(summary(lm(cl1Y_Fact ~ cl1X1_Fact + cl1X2_Fact - 1))$sigma)
  
  #premiere composante de G a normer par le residu de la regression pour gagner du temps
  regY_Fact <-lm(Y - T%*%ini_D ~ cl1Y_Fact - 1 )
  ini_B_Fact = matrix(coef(regY_Fact), nrow=1,ncol = qY, byrow=TRUE)
  sig_lm_qY_Fact <- rep(NA,qY)
  for(i in 1:qY){
    sig_lm_qY_Fact[i] <- summary(regY_Fact)[[i]]$sigma
  }
  ini_sigma2Y_Fact = mean(sig_lm_qY_Fact^2) #initialisation du parametre sigma2Y 
  #Regression de l'initialisation de F sur celle de G pour estimation C_sim 
  regG_Fact <-lm(cl1Y_Fact ~ cl1X1_Fact + cl1X2_Fact -1)
  ini_C1_Fact <- regG_Fact$coefficients[1] #initialisation du parametre C1
  ini_C2_Fact <- regG_Fact$coefficients[2] #initialisation du parametre C2
  
  
  ############################################
  # Debut initialisation
  ############################################
  A1 = ini_A1_Fact
  A2 = ini_A2_Fact
  B = ini_B_Fact
  C1 = ini_C1_Fact 
  C2 = ini_C2_Fact 
  D = ini_D #New
  D1 = ini_D1 #New
  D2 = ini_D2 #New
  sigma2_X1_chap = ini_sigma2X1_Fact 
  sigma2_X2_chap = ini_sigma2X2_Fact 
  sigma2_Y_chap = ini_sigma2Y_Fact
  #D'o?
  Psi_X1 = diag(sigma2_X1_chap, qX1)#necessaire pour ds boucle calculer M-i et Sigma_i de la loi qui donne F_tilde, G_tilde,Phi_tilde et Gamma_tilde
  Psi_X2 = diag(sigma2_X2_chap, qX2)
  Psi_Y = diag(sigma2_Y_chap, qY)
  ########################
  # 
  # FIN INITIALISATION            
  #            
  ########################
  
  
  
  
  
  
  it = 0 # ATTENTION it = 0 cree des objets de taille nulle et des bugs ds while{} a partir ligne 118 si it ne passe pas a it=0+1
  
  #On cree un vecteur qui va contenir les differences d'une iteration a l'autre pour voir l'evolution de la convergence
  diffgraph =  matrix(rep(0,(nb_it-1)), ncol=1)
  
  diff = 1 # valeur initiale du parametre mesurant le changemet de valeur de theta^[t] d'une iteration a l'autre
  #D'ou, l'initialisation se traduit au niveau de diffgraph par : 
  diffgraph[1]=1
  
  
  detE2 =  matrix(rep(0,nb_it), ncol=1)
  
  while( diff > epsilon && (it = it+1) < nb_it)#en utilisant &&, seule la premiere partie de la condition est evaluee
    #et quand elle est vraie it passe a it+1 ET l'instruction est executee
  {
    #Calcule pour les n individus i de : phi_tilde, gama_tilde, f_tilde et g_tilde :
    #Besoin d'abord de definir M_i et GAMMA_i resp. Le vecteur moy et la matrice de var-cov de la loi normale de H_i|Z_i 
    ###################################################################################################################
    ##                 Def d'elements necessaires a la def des param a l'iteration [t+1] :                           ##
    ###################################################################################################################  
    
    
    #Quand Les param sont initialises sous forme matrix est que dim theorique conservees : on utilise ce code de E1 et E2
    #ou codes comme dans theorie sans les t() utilises quand on les initialise en numeric pour recuperer les dim theoriques
    E1_1 = cbind( ((C1)^2 + (C2)^2 + 1)*B
                  , C1*A1
                  ,C2*A2# c1 et c2 reels dans le cas de un facteur d'ou "*" et non "%*%"
    )
    E1_2 = cbind( C1*B
                  ,A1
                  ,t(rep(0,qX2))
    )
    E1_3 = cbind( C2*B                
                  ,t(rep(0,qX1))
                  ,A2
    )
    
    E1 = rbind(E1_1, E1_2, E1_3)
    
    
    E2_1 = cbind( ((C1)^2 + (C2)^2 + 1)*t(B)%*%B + Psi_Y
                  , C1*t(B)%*%A1
                  ,C2*t(B)%*%A2# c1 et c2 r?els dans le cas de un facteur d'o? "*" et non "%*%"
    )
    E2_2 = cbind( C1*t(A1)%*%B
                  ,t(A1)%*%A1 + Psi_X1
                  ,matrix(rep(0,qX1*qX2),ncol=qX2)
    )
    E2_3 = cbind( C2*t(A2)%*%B                
                  ,matrix(rep(0,qX1*qX2),ncol=qX1)
                  ,t(A2)%*%A2 + Psi_X2
    )
    
    E2 = rbind(E2_1, E2_2, E2_3) #superpose les lignes E_1 E_2 puis E_3 pour construire la matrice E ie : D2
    
    detE2[it,] = det(E2)
    
    
    E3 = array( 1:((qY+qX1+qX2)*1*n), dim = c(qY+qX1+qX2, 1, n))
    #New
    for(i in 1:n){
      E3[,,i] = c(   if(qT==1){Y[i,] - t(D) * T[i,]}else{ Y[i,] - t(D) %*% T[i,]},#car D, D1 et D2 sont des matrices
                     if(qT1==1){X1[i,] - t(D1) * T1[i,]}else{ X1[i,] - t(D1) %*% T1[i,]},
                     if(qT2==1){X2[i,] - t(D2) * T2[i,]}else{ X2[i,] - t(D2) %*% T2[i,]}
      )
    }
    
    #M_i correspond a E1%*%solve(E2)%*%E3[,,i]
    
    #Codons GaMMA_i
    
    E4 = matrix(
      c( (C1)^2 + (C2)^2 + 1, 
         C1,
         C2,
         C1,
         1,
         0,
         C2,
         0,
         1
      )
      , ncol=kG+kF1+kF2)
    
    
    GAMMA = E4 - E1%*%ginv(E2)%*%t(E1)  
    
    
    
    
    EE = E1%*%ginv(E2)
    tEE = t(ginv(E2))%*%t(E1)
    fg = prod_tab_mat(EE , E3) #fg est un tableau de i=1 a n elements fg[,,i]
    #chaque fg[,,i] contient, les kG premieres ligne m_{1i} c'est a dire le premiere element moyenne de h_i|z_i
    #et les kF1 lignes suivantes contient m_{2i} et les kF2 dernieres lignes contient m_{3i}
    
    
    #G_tilde
    G_tilde = array( fg[1:kG,,] , dim = c(kG , dim(fg)[2] , dim(fg)[3]) )
    #G_tilde_bar
    G_tilde_bar = apply(G_tilde, 1:2, mean)
    
    #F1_tilde
    F1_tilde = array( fg[(kG+1):(kG+kF1),,] , dim = c(kF1 , dim(fg)[2] , dim(fg)[3]) )
    #F1_tilde_bar
    F1_tilde_bar = apply(F1_tilde, 1:2, mean)
    
    #F2_tilde
    F2_tilde = array( fg[(kG+kF1+1):(kG+kF1+kF2),,] , dim = c(kF2 , dim(fg)[2] , dim(fg)[3]) )
    #F2_tilde_bar
    F2_tilde_bar = apply(F2_tilde, 1:2, mean)
    
    
    
    #G_tilde_Y_bar
    G_tilde_Y_bar <- t(moy_TildeData(G_tilde, Y))
    
    #G_tilde_T_prim_bar
    G_tilde_T_prim_bar <- moy_TildeData(G_tilde, T ) #New
    #Rq : on ne met pas t(T) car code ds la fction moy_TildeData()
    
    #G_tilde_T_bar
    G_tilde_T_bar <- t(moy_TildeData(G_tilde, T )) #New
    
    
    
    
    #F1_tilde_X1_bar
    F1_tilde_X1_bar <- t(moy_TildeData(F1_tilde, X1))
    
    #F2_tilde_X2_bar
    F2_tilde_X2_bar <- t(moy_TildeData(F2_tilde, X2))
    
    #F1_tilde_T1_prim_bar
    F1_tilde_T1_prim_bar <- moy_TildeData(F1_tilde, T1)
    
    #F2_tilde_T2_prim_bar
    F2_tilde_T2_prim_bar <- moy_TildeData(F2_tilde, T2)
    
    #F1_tilde_T1_bar
    F1_tilde_T1_bar <- t(moy_TildeData(F1_tilde, T1))
    
    #F2_tilde_T2_bar
    F2_tilde_T2_bar <- t(moy_TildeData(F2_tilde, T2))
    
    #Y_T_prim_bar
    Y_T_prim_bar <- moy_T_Prim_Data(Y, T)
    
    #X1_T1_prim_bar
    X1_T1_prim_bar <- moy_T_Prim_Data(X1, T1)
    
    #X2_T2_prim_bar
    X2_T2_prim_bar <- moy_T_Prim_Data(X2, T2)
    
    
    
    
    moy_G_tilde = apply(G_tilde, 1:2, mean)#objet cree pour creer l'objet G_tilde_bar_2 qui va suivre
    
    moy_F1_tilde = apply(F1_tilde, 1:2, mean)#objet cree pour creer l'objet F1_tilde_bar_2 qui va suivre
    
    moy_F2_tilde = apply(F2_tilde, 1:2, mean)#objet cree pour creer l'objet F2_tilde_bar_2 qui va suivre
    
    
    #G_tilde_bar_2
    G_tilde_bar_2 <- moy_G_tilde%*%t(moy_G_tilde)
    
    #F1_tilde_bar_2
    F1_tilde_bar_2 <- moy_F1_tilde%*%t(moy_F1_tilde)
    
    #F2_tilde_bar_2
    F2_tilde_bar_2 <- moy_F2_tilde%*%t(moy_F2_tilde)
    
    
    
    #F1G_tilde_bar
    F1G_tilde_bar <- moy_fg_tilde( F1_tilde , G_tilde )
    
    #F2G_tilde_bar
    F2G_tilde_bar <- moy_fg_tilde( F2_tilde , G_tilde )
    
    #F1F2_tilde_bar
    F1F2_tilde_bar <- moy_fg_tilde( F1_tilde , F2_tilde )
    
    
    #GAMMA_12_bar
    GAMMA_12_bar <- GAMMA[1:kG, (kG+1):(kG+kF1)]
    
    #GAMMA_13_bar
    GAMMA_13_bar<- GAMMA[1:kG, (kG+kF1+1):(kG+kF1+kF2)]
    
    #GAMMA_23_bar
    GAMMA_23_bar <- GAMMA[(kG+1):(kG+kF1), (kG+kF1+1):(kG+kF1+kF2)]
    
    
    #Gamma_tilde
    Gamma_tilde = array( rep(0) , dim = c( kG , kG , n))
    for(i in 1:n){
      Gamma_tilde[,,i] = ( (EE%*%E3[,,i])^2 )[1:kG,] + GAMMA[1:kG, 1:kG]
    }
    
    #Phi1_tilde
    Phi1_tilde = array( rep(0) , dim = c( kF1 , kF1 , n))
    for(i in 1:n){
      Phi1_tilde[,,i] = ( (EE%*%E3[,,i])^2 )[(kG+1):(kG+kF1),] + GAMMA[(kG+1):(kG+kF1), (kG+1):(kG+kF1)]
    }
    
    #Phi2_tilde
    Phi2_tilde = array( rep(0) , dim = c( kF2 , kF2 , n))
    for(i in 1:n){
      Phi2_tilde[,,i] = ( (EE%*%E3[,,i])^2 )[(kG+kF1+1):(kG+kF1+kF2),] + GAMMA[(kG+kF1+1):(kG+kF1+kF2), (kG+kF1+1):(kG+kF1+kF2)]
    }
    
    
    #Gamma_tilde_bar
    Gamma_tilde_bar <- apply(Gamma_tilde, 1:2, mean)
    
    #Phi1_tilde_bar
    Phi1_tilde_bar <- apply(Phi1_tilde, 1:2, mean)
    
    #Phi2_tilde_bar
    Phi2_tilde_bar <- apply(Phi2_tilde, 1:2, mean)
    
    #Phi1Phi2_tilde_bar
    Phi1Phi2_tilde_bar <- moy_fg_tilde( Phi1_tilde , Phi2_tilde )
   
    
    #T_T_prim_bar :
    T_T_prim_bar = moy_T_T_Prim(T,T)
    
    #T1_T1_prim_bar :
    T1_T1_prim_bar = moy_T_T_Prim(T1,T1)
    
    #T2_T2_prim_bar :
    T2_T2_prim_bar = moy_T_T_Prim(T2,T2)
    
    
    
    ###################################################################################################################
    ##               Fin def d'elements necessaires a la def des param a l'iteration [t+1]                           ##
    ###################################################################################################################
    
    ######################################################################################################################
    ###     A ACTUALISER SOUS FM3 
    
    #Stock des valeurs des parametres a l'iteration it avant actualisation et l'obtention de leur valeur pour it+1
    #Objectif = nous permettre de calculer la difference de valeur d'une iteration a l'autre pour le critere d'arret
    #THETA[it] :
    A1_it = A1
    A2_it = A2
    B_it = B
    C1_it = C1
    C2_it = C2
    D_it = D #New
    D1_it = D1 #New
    D2_it = D2 #New
    sigma2_X1_chap_it =  sigma2_X1_chap
    sigma2_X2_chap_it =  sigma2_X2_chap
    sigma2_Y_chap_it = sigma2_Y_chap
    
    ##############################################################
    #THETA[it+1,] 
    #Parametres a l'iteration it+1 :
    B = t((G_tilde_Y_bar - Y_T_prim_bar%*%solve(T_T_prim_bar)%*%G_tilde_T_bar )%*%
            solve(Gamma_tilde_bar - G_tilde_T_prim_bar%*%solve(T_T_prim_bar)%*%G_tilde_T_bar))
    
    A1 = t((F1_tilde_X1_bar - X1_T1_prim_bar%*%solve(T1_T1_prim_bar)%*%F1_tilde_T1_bar )%*%
             solve(Phi1_tilde_bar - F1_tilde_T1_prim_bar%*%solve(T1_T1_prim_bar)%*%F1_tilde_T1_bar))
    
    A2 = t((F2_tilde_X2_bar - X2_T2_prim_bar%*%solve(T2_T2_prim_bar)%*%F2_tilde_T2_bar )%*%
             solve(Phi2_tilde_bar - F2_tilde_T2_prim_bar%*%solve(T2_T2_prim_bar)%*%F2_tilde_T2_bar))
    
    C1 = solve( Phi1Phi2_tilde_bar - (GAMMA_23_bar + F1F2_tilde_bar)^2  #denominateur   
    )%*%((GAMMA_12_bar + F1G_tilde_bar)*Phi2_tilde_bar - (GAMMA_13_bar + F2G_tilde_bar)*(GAMMA_23_bar + F1F2_tilde_bar))#numerateur
    C2 = solve( Phi1Phi2_tilde_bar - (GAMMA_23_bar + F1F2_tilde_bar)^2   #denominateur
    )%*%((GAMMA_13_bar + F2G_tilde_bar)*Phi1_tilde_bar - (GAMMA_12_bar + F1G_tilde_bar)*(GAMMA_23_bar + F1F2_tilde_bar))#numerateur
    
    D = t((Y_T_prim_bar - t(B)%*%G_tilde_T_prim_bar)%*%solve(T_T_prim_bar))  
    D1 = t((X1_T1_prim_bar - t(A1)%*%F1_tilde_T1_prim_bar)%*%solve(T1_T1_prim_bar))
    D2 = t((X2_T2_prim_bar - t(A2)%*%F2_tilde_T2_prim_bar)%*%solve(T2_T2_prim_bar))
    
    ####################################################
    #    sigma2_Y_chap
    #######################################################
    for(i in 1:n){
      dud1_Y[i,] <- Y[i,]-t(D)%*%T[i,]
    }
    res1_Y <- apply((dud1_Y)^2,1,sum)
    res2_Y <- c(sum(B^2)*Gamma_tilde)
    tab_res3_Y = array(rep(0), dim=c(n,1,n))
    for(i in 1:n){
      tab_res3_Y[,,i] = dud1_Y[i,]%*%t(B)%*%G_tilde[,,i]
    }
    res3_Y <- apply(tab_res3_Y, 1:2, mean)
    res_Y_chap <- res1_Y + res2_Y - 2*res3_Y
    sigma2_Y_chap <- max((1/(n*qY))*sum(res_Y_chap))
    ####################################################
    #    sigma2_X1_chap
    ####################################################
    for(i in 1:n){
      dud1_X1[i,] <- X1[i,]-t(D1)%*%T1[i,]
    }
    res1_X1 <- apply((dud1_X1)^2,1,sum)
    res2_X1 <- c(sum(A1^2)*Phi1_tilde)
    tab_res3_X1 = array(rep(0), dim=c(n,1,n))
    for(i in 1:n){
      tab_res3_X1[,,i] = dud1_X1[i,]%*%t(A1)%*%F1_tilde[,,i]
    }
    res3_X1 <- apply(tab_res3_X1, 1:2, mean)
    #Pour un idividu i, ca revient au meme que  ce qui semble plus logique : t((X[i,]-mu_X[,,it]))%*%(as.matrix(A[,,it]))*F_tilde[,,i]
    res_X1_chap <- res1_X1 + res2_X1 - 2*res3_X1
    sigma2_X1_chap<-max((1/(n*qX1))*sum(res_X1_chap))
    ####################################################
    #    sigma2_X2_chap
    ####################################################
    for(i in 1:n){
      dud1_X2[i,] <- X2[i,]-t(D2)%*%T2[i,]
    }
    res1_X2 <- apply((dud1_X2)^2,1,sum)
    res2_X2 <- c(sum(A2^2)*Phi2_tilde)
    tab_res3_X2 = array(rep(0), dim=c(n,1,n))
    for(i in 1:n){
      tab_res3_X2[,,i] = dud1_X2[i,]%*%t(A2)%*%F2_tilde[,,i]
    }
    res3_X2 <- apply(tab_res3_X2, 1:2, mean)
    #Pour un idividu i, ca revient au meme que  ce qui semble plus logique : t((X[i,]-mu_X[,,it]))%*%(as.matrix(A[,,it]))*F_tilde[,,i]
    res_X2_chap <- res1_X2 + res2_X2 - 2*res3_X2
    sigma2_X2_chap<-max((1/(n*qX2))*sum(res_X2_chap))
    
    
    
    ########################
    #  Psi_X1, Psi_X2 et Psi_Y a it+1
    ########################
    Psi_Y <- sigma2_Y_chap*diag(1,qY)
    Psi_X1 <- sigma2_X1_chap*diag(1,qX1)
    Psi_X2 <- sigma2_X2_chap*diag(1,qX2)
    
    
    ##############################
    #---- Conversion des parametres de matrix a numeric----
    #On convertit C1 et C2 en numeric car apres estimation ce sont des matrix et le debut du code est cree avec des C1 et C2 numeric.
    
    C1 = as.numeric(C1)
    C2 = as.numeric(C2)
    
    #car a l'initialisation ceux sont des numeric et apres application de la boucle while il deviennent des
    #matrix. Si on ne le fait pas, l'application de la boucle lors de la deuxieme fois ne peut etre correcte.
    #Rq : elle n'est pas appliquee aux sigma2 car apres appli de la boucle while leur class n'est pas changee
    
    
    
    
    #mesure du changement d'iteration :
    diff = sum(
      sum(abs(D -  D_it)/abs(D)),
      sum(abs(D1 -  D1_it)/abs(D1)),
      sum(abs(D2 -  D2_it)/abs(D2)),
      sum(abs( B -  B_it)/abs(B)) ,
      sum(abs( A1 -  A1_it)/abs(A1)),
      sum(abs( A2 -  A2_it)/abs(A2)) ,
      (abs( C1 -  C1_it)/abs(C1)) ,
      (abs( C2 -  C2_it)/abs(C2)) ,
      (abs(sigma2_Y_chap - sigma2_Y_chap_it)/abs(sigma2_Y_chap)) ,
      (abs(sigma2_X1_chap - sigma2_X1_chap_it)/abs(sigma2_X1_chap)) ,
      (abs(sigma2_X2_chap - sigma2_X2_chap_it)/abs(sigma2_X2_chap))
    )
    
    diffgraph[it+1]=diff
    
    
    cat('iteration', it, '\n', '\n',
        'difference', diffgraph[it], 
        '\n','D_chap', D,
        '\n','D1_chap', D1,
        '\n','D2_chap', D2,
        '\n','B_chap', B,
        '\n','A1_chap', A1, '\n','A2_chap', A2,
        '\n','C1_chap', C1, '\n','C2_chap', C2,
        '\n', 'sigma2_Y_chap', sigma2_Y_chap,
        '\n', 'sigma2_X1_chap', sigma2_X1_chap,
        '\n', 'sigma2_X2_chap', sigma2_X2_chap,  '\n'
    )
    
  }#fin boucle while du code
  
  
  return(list(Factors=cbind(G_tilde,F1_tilde,F2_tilde),
              C=c(C1,C2),
              B= B,
              A1= A1,
              A2= A2,
              D=D,
              D1=D1,
              D2=D2,
              sigma2 = c(sigma2_Y_chap, sigma2_X1_chap, sigma2_X2_chap),
              Diff = diff,
              Diff_it=diffgraph,
              Conv_it=it
  )
  )
  
  
}#end function SEM_FM3
