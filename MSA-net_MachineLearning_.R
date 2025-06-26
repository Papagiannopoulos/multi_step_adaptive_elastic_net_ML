
# To run this script, the user must provide three input files with the following names:
#> 1) y
#> 2) metaboRank
#> 3) abr_adiposity
#> 
# Libraries
library(readr); library(readxl); library(xlsx)
library(caret); library(msaenet)#MULTI STEP ADAPTIVE ELASTIC NET
library(dplyr); library(purrr)
library(ggplot2); library(ggpubr); 
library(doParallel); library(parallel); library(doRNG)#For Parallelization
library(MASS); library(splitstackshape)

#metaboRank <- read_excel("")
#y <- read_excel("")

# 1) y: Adiposity Data Frame (Y: Feature)
#str(y)
# data.frame:	70591 obs. of  8 variables:
# Body Fat          : num  0.983 0.137 1.073 0.317 -0.133 ...
# Waist             : num  -0.8231 -0.0584 0.5151 0.993 -0.7275 ...
# Hip               : num  -1.876 0.294 0.873 1.018 0.149 ...
# Wasit to Hip Ratio: num  0.6341 -0.3614 -0.0101 0.5854 -1.3039 ...
# BMI               : num  -0.619 -0.431 0.296 0.228 -0.626 ...
# ABSI              : num  -0.00178 1.0246 0.41346 1.41344 -0.17265 ...
# HI                : num  -1.672 1.958 0.981 1.169 1.915 ...
# WHI               : num  1.206 -0.196 -0.229 0.551 -1.263 ...
# .
# .
# .

# Structure of the used data sets
# 2) metaboRank: Metabolomics Data Frame (X: Features)
#str(metaboRank)
# data.frame:	~150,000 obs. of  250 variables:
# my_var    : num  0.983 0.137 1.073 0.317 -0.133 ...
# X23400.0.0: num  1.3303 0.0362 -0.0591 -0.8152 1.4455 ...
# X23401.0.0: num  0.93 -0.111 -0.342 -0.892 1.363 ...
# X23402.0.0: num  1.12 -0.245 -0.432 -0.78 1.121 ...
# X23403.0.0: num  0.689 -0.657 -0.571 -0.602 0.773 ...
# X23404.0.0: num  0.8893 0.0822 -0.3364 -1.0878 1.513 ...
# X23405.0.0: num  0.7212 0.0154 -0.2438 -0.9506 1.5361 ...
# X23406.0.0: num  1.538 0.615 0.916 0.172 0.729 ...
# X23407.0.0: num  -0.0708 -1.378 -0.5791 0.2737 0.0938 ...
# .
# .
# .

# 3) Abbreviations
# abr_adiposity <- c("BF%","WC","HC","WHR","BMI","ABSI","HI", "WHI")

# For clarity, simulated data is provided to demonstrate the required input formats.
# Replace the simulated data with your own data in the specified format before running the code.

set.seed(4423)
mu <- rep(0, 4)
Sigma <- var(matrix(c(1, .1, .9, .9, .1, 1, .9, .9, .9, .9, 1, .9, .9, .9, .9, 1), 4, 4))
metaboRank <- data.frame(mvrnorm(10^3, mu, Sigma))

y <- metaboRank[,c(1,1)] + rnorm(10^3, 1, 1)

metaboRank <- data.frame(my_var = y[,1], metaboRank)

abr_adiposity <- paste0(1:length(y), "X")

##### Stability Selection #####
#> If the code crashes at this step, the user may need to adjust the hyperparameter ranges.
#> 
set.seed(3487)
# Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 23423)

features <- list()
N <- 100 # Number of Splits
# Hyperparameters
steps <- 5
criteria <- c("lambda.1se","lambda.min")
gamma <- c(0.1,0.5,1,1.5,2,2.5,3)
for(p in 1:length(y)){
  print(p)
  set.seed(572920)
  features[[p]] <- vector()
  metaboRank$my_var <- y[,p]
  for(j in 1:N){
    print(j)
    train_rows <- metaboRank$my_var %>% createDataPartition(p = 0.80, list = F)
    train <- metaboRank[train_rows,]
    
    x1 <- sample(criteria, 1, replace = T)
    x2 <- sample(gamma, 1, replace = T)
    print(x1); print(x2)
    msanet <- msaenet(as.matrix(train[,-match("my_var",names(train))]),train[,match("my_var",names(train))],
                      family = "gaussian",init = "enet",tune = "cv",nfolds = 10,
                      rule = x1,alphas = seq(0.1,1,0.1),verbose = T,
                      parallel = T,scale = x2,seed = T,nsteps = steps,tune.nsteps = "ebic")
    
    msanet_coef <- data.frame(names = rownames(msanet$beta),coef1 = msanet$beta[,1])
    msanet_coef <-  msanet_coef[msanet_coef$coef1 != 0,]
    features[[p]] <- c(features[[p]], msanet_coef$names)
  }
  features[[p]] <- data.frame(Feature = names(table(features[[p]])), Counts = as.numeric(table(features[[p]])))
  features[[p]] <- arrange(features[[p]], desc(Counts))

  if(!dir.exists("FS")){ dir.create("FS")}
  write.xlsx(features[[p]], paste0("FS/",names(y)[p],"_Feature_Selection_for.xlsx"), row.names = F, showNA = F)
}
stopCluster(cl)

##### Tuning #####
set.seed(34587)
# Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 23423)

msanet <- list()
msanet_rmse <- list()
msanet_pearson <- list()
test <- list(); tune <- list()
for(p in 1:length(y)){
  set.seed(17997)
  tune[[p]] <- matrix(NA,ncol = length(criteria), nrow = length(gamma))
  colnames(tune[[p]]) <- criteria
  rownames(tune[[p]]) <- gamma
  names(tune)[p] <- names(y)[p]
  
  features <- read.xlsx(paste0("FS/",names(y)[p],"_Feature_Selection_for.xlsx"), sheetIndex = 1)
  for(c in 1:length(criteria)){
    for(g in 1:length(gamma)){
      print(paste0(names(y)[p],", criteria = ",criteria[c],", gamma = ",gamma[g]))
      
      metaboRank$my_var <- y[,p]
      train_rows <- metaboRank$my_var %>% createDataPartition(p = 0.80, list = F)
      train <- metaboRank[train_rows,c(1,match(features$Feature[features$Counts >= N*0.6], names(metaboRank)))]
      test <- metaboRank[-train_rows,c(1,match(features$Feature[features$Counts >= N*0.6], names(metaboRank)))]
      
      msanet <- msaenet(as.matrix(train[,-match("my_var",names(train))]),train[,match("my_var",names(train))],
                        family = "gaussian",init = "enet",tune = "cv",nfolds = 10,
                        rule = criteria[c],alphas = seq(0.1,1,0.1),verbose = T,
                        parallel = T,scale = gamma[g],seed = T,nsteps = steps,tune.nsteps = "max")
      
      msanet_rmse <- RMSE(test[,match("my_var",names(test))],predict(msanet,newx=as.matrix(test[,-match("my_var",names(test))])))
      msanet_pearson <- cor.test(test[,match("my_var",names(test))],predict(msanet,newx=as.matrix(test[,-match("my_var",names(test))])))$est
      
      tune[[p]][g,c] <- paste0(round(msanet_rmse,3),"/",round(msanet_pearson,3),"/",msanet$model$df,"/",round(msanet_rmse/msanet$model$df,4))
    }
  }
}
stopCluster(cl)

# Save Tuning results
save_tuning <- data.frame()
tune_parametrs <- list()
for(p in 1:length(tune)){
  d <- as.data.frame(tune[[p]])
  d1 <- as.matrix(cSplit(d,names(d), sep = "/", type.convert = F))
  colnames(d1) <- c("lsd_rmse","lsd_pearson","lsd_df","lsd_ratio","lmin_rmse","lmin_pearson","lmin_df","lmin_ratio")
  rownames(d1) <- gamma
  d1 <- apply(d1,2,as.numeric)
  
  d1[,4] <- d1[,1] / d1[,3]
  d1[,8] <- d1[,5] / d1[,7]
  print(d1)
  tune_parametrs[[p]] <- which(d1[,grepl("ratio", colnames(d1))] == min(d1[,grepl("ratio", colnames(d1))]),arr.ind = T)
  
  save_tuning <- rbind(save_tuning,data.frame(d1,var = names(y)[p]))
}
write.xlsx(save_tuning, "full_tuning.xlsx",row.names = F, showNA = F, col.names = T)

##### Tuning - Figures #####
save_tuning <- as.data.frame(read_excel("full_tuning.xlsx"))
save_tuning <- reshape(save_tuning, direction = "long", idvar = "id", timevar = "criteria", varying = list(c(1,5),c(2,6),c(3,7),c(4,8)),
                       v.names = c("rmse","pearson","df","ratio"), sep = "_")
save_tuning <- save_tuning[!grepl("Resid",save_tuning$var),]
save_tuning <- data.frame(save_tuning,gamma = rep(gamma,length(y)))
save_tuning$criteria <- factor(save_tuning$criteria, c(1,2),c("Lambda 1sd","Lambda min"))
# Candidate Hyperparameters
# After reviewing the tuning visualizations, the user should select the optimal hyperparameters
# based on RMSE and Pearson metrics, assign their positions in the "row" and "col" variables,
# then rerun the script to generate and save the updated graphs.
row <- sample(1:length(gamma), length(y), replace = T) # refer to gamma hyperparameters position. Its length equals the length of Y feature.
col <-  sample(1:length(criteria), length(y), replace = T) # refer to criteria hyperparameters position. Its length equals the length of Y feature.

gtune <- list()
for(p in 1:length(y)){
  g1 <- ggplot(save_tuning[save_tuning$var %in% names(y)[p],], aes(x=gamma, y=rmse, colour = criteria)) + geom_point(cex = 5) + 
    labs(x="", y="RMSE", title = abr_adiposity[p])+ 
    theme_bw()+ theme(plot.title = element_text(size=15),axis.text.x= element_text(size=8),
                      axis.text.y= element_text(size=12), axis.title=element_text(size=12))+
    geom_line(lwd = 1.5)+scale_x_continuous(breaks = gamma, limits = c(0,3.3))+geom_text(aes(label = df), size = 5, col = "red")+
    guides(size = "none", color = "none")+ylim(c(min(save_tuning[save_tuning$var %in% names(y)[p],]$rmse,na.rm = T)*0.95,
                                                 max(save_tuning[save_tuning$var %in% names(y)[p],]$rmse,na.rm = T)*1.02))+
    annotate("label",x = save_tuning[save_tuning$var %in% names(y)[p],][save_tuning[save_tuning$var %in% names(y)[p],]$criteria %in% levels(save_tuning$criteria)[col[p]],8][row[p]],
             y= save_tuning[save_tuning$var %in% names(y)[p],][save_tuning[save_tuning$var %in% names(y)[p],]$criteria %in% levels(save_tuning$criteria)[col[p]],3][row[p]],
             label = "X", size = 5, col = c("red1","cyan4")[col[p]], vjust = 2.7)
  
  g2 <- ggplot(save_tuning[save_tuning$var %in% names(y)[p],], aes(x=gamma, y=pearson, colour = criteria)) + geom_point(cex = 5) + 
    labs(x="", y="Pearson",title = abr_adiposity[p])+ 
    theme_bw()+ theme(plot.title = element_text(size=15),axis.text.x= element_text(size=8),
                      axis.text.y= element_text(size=12), axis.title=element_text(size=12))+
    geom_line(lwd = 1.5)+scale_x_continuous(breaks = gamma, limits = c(0,3.3))+guides(size = "none", color = "none")+
    geom_text(aes(label = df, size = 5), col = "red")+
    ylim(c(min(save_tuning[save_tuning$var %in% names(y)[p],]$pearson)*0.97,
           max(save_tuning[save_tuning$var %in% names(y)[p],]$pearson)*1.03))+
    annotate("label",x = save_tuning[save_tuning$var %in% names(y)[p],][save_tuning[save_tuning$var %in% names(y)[p],]$criteria %in% levels(save_tuning$criteria)[col[p]],8][row[p]],
             y= save_tuning[save_tuning$var %in% names(y)[p],][save_tuning[save_tuning$var %in% names(y)[p],]$criteria %in% levels(save_tuning$criteria)[col[p]],4][row[p]],
             label = "X",size = 5, col = c("red1","cyan4")[col[p]], vjust=2.7)
  
  if(p %in% c(4,8)){
    gtune[[p]] <- ggarrange(g1+xlab("Hyperparameter Gamma"),g2+xlab("Hyperparameter Gamma"))
  }else{
    gtune[[p]] <- ggarrange(g1,g2)
  }
  gtune[[p]]
}
#> The length of the following ggarrange plot should equals the length of y element
ggarrange(plotlist = gtune)
ggsave("hypeparametres.png", dpi = 1000, width = 12, height = 8)

##### Final Model & Save Results #####
set.seed(67489)
# Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 9401)

msanet <- list()
msanet_rmse <- list()
msanet_pearson <- list()
test <- list()
metrics <- list()
N <- 100 # Bootstrapping
coeff <- list()
for(p in 1:length(y)){
  print(p)
  msanet_rmse[[p]] <- NA
  msanet_pearson[[p]] <- NA
  coeff[[p]] <- data.frame()
  
  features <- read.xlsx(paste0("FS/",names(y)[p],"_Feature_Selection_for.xlsx"), sheetIndex = 1)
  metaboRank$my_var <- y[,p]
  for(i in 1:N){
    print(paste0("N = ", i))
    train_rows <- metaboRank$my_var %>% createDataPartition(p = 0.80, list = F)
    train <- metaboRank[train_rows,c(1,match(features$Feature[features$Counts >= N*0.6], names(metaboRank)))]
    test <- metaboRank[-train_rows,c(1,match(features$Feature[features$Counts >= N*0.6], names(metaboRank)))]
    
    msanet[[p]] <- msaenet(as.matrix(train[,-match("my_var",names(train))]),train[,match("my_var",names(train))],
                           family = "gaussian",init = "enet",tune = "cv",nfolds = 10,
                           rule = criteria[col[p]],alphas = seq(0.1,1,0.1),verbose = T,
                           parallel = T,scale = gamma[row[p]],seed = T,nsteps = steps,tune.nsteps = "max")
    
    # Metrics
    msanet_rmse[[p]] <- c(msanet_rmse[[p]], RMSE(test[,match("my_var",names(test))],predict(msanet[[p]],newx=as.matrix(test[,-match("my_var",names(test))]))))
    msanet_pearson[[p]] <- c(msanet_pearson[[p]],cor.test(test[,match("my_var",names(test))],predict(msanet[[p]],newx=as.matrix(test[,-match("my_var",names(test))])))$est)
    # Coefficients
    cc <- data.frame(names = rownames(msanet[[p]]$beta),coef1 = msanet[[p]]$beta[,1])
    cc <-  cc[cc$coef1 != 0,]
    
   if(i == 1){
     coeff[[p]]  <- data.frame(cc)
     names(coeff[[p]])[length(coeff[[p]])] <- paste0("Iter",i)
   }else{
     coeff[[p]]  <- merge(coeff[[p]], cc, by = "names", all = T)
     names(coeff[[p]])[length(coeff[[p]])] <- paste0("Iter",i)
   }
}

  if(!dir.exists("coef_profiles")){ dir.create("coef_profiles")}
  write.xlsx(coeff[[p]], paste0("coef_profiles/Coef_N100_",names(y)[p],".xlsx"), showNA = F)
  write.xlsx(msanet_rmse[[p]], paste0("coef_profiles/RMSE_N100_",names(y)[p],".xlsx"), showNA = F)
  write.xlsx(msanet_pearson[[p]], paste0("coef_profiles/Pearson_N100_",names(y)[p],".xlsx"), showNA = F)
}
stopCluster(cl)
