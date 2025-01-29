
##### Stability Selection #####
set.seed(3487)
#Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 23423)

features <- list()
N <- 100#Number of Splits
#Hyperparameters
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
  features[[p]]$Abr <- UKB_abr$Metabolite[match(features[[p]]$Feature, UKB_abr$data_field)]
  
  write.xlsx(features[[p]], paste0("FS/",names(y)[p],"_Feature_Selection_for.xlsx"), row.names = F, showNA = F)
}
stopCluster(cl)

##### Tuning #####
set.seed(34587)
#Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 23423)

msanet <- list()
msanet_rmse <- list()
msanet_pearson <- list()
test <- list(); tune <- list()
#Hyperparameters
steps <- 5
criteria <- c("lambda.1se","lambda.min")
gamma <- c(0.1, 0.5, 1, 1.5, 1.8, 2, 2.5, 3)
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
      train <- metaboRank[train_rows,c(1,match(features$Feature[features$Counts >= 60], names(metaboRank)))]
      test <- metaboRank[-train_rows,c(1,match(features$Feature[features$Counts >= 60], names(metaboRank)))]
      
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

#Save Tuning results
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
write.xlsx(save_tuning, "full_tuning_male.xlsx",row.names = F, showNA = F, col.names = T)

##### Tuning - Figures #####
save_tuning <- as.data.frame(read_excel("full_tuning_male.xlsx"))
save_tuning <- reshape(save_tuning, direction = "long", idvar = "id", timevar = "criteria", varying = list(c(1,5),c(2,6),c(3,7),c(4,8)),
                       v.names = c("rmse","pearson","df","ratio"), sep = "_")
save_tuning <- save_tuning[!grepl("Resid",save_tuning$var),]
save_tuning <- data.frame(save_tuning,gamma = rep(gamma,length(y)))
save_tuning$criteria <- factor(save_tuning$criteria, c(1,2),c("Lambda 1sd","Lambda min"))
#Candidate Hyperparameters
row <- c(6,4,7,7,5,6,8,6)
col <- c(1,1,2,1,1,2,2,2)
gtune <- list()
for(p in 1:length(y)){
  g1 <- ggplot(save_tuning[save_tuning$var %in% names(y)[p],], aes(x=gamma, y=rmse, colour = criteria)) + geom_point(cex = 5) + 
    labs(x="", y="RMSE", title = abr_adiposity[p])+ #guides(fill = "", colour = "")+
    theme_bw()+ theme(plot.title = element_text(size=15),axis.text.x= element_text(size=8),
                      axis.text.y= element_text(size=12), axis.title=element_text(size=12))+
    geom_line(lwd = 1.5)+scale_x_continuous(breaks = gamma, limits = c(0,3.3))+geom_text(aes(label = df), size = 5, col = "red")+#, hjust = rep(c(-2,2),each=7), vjust = rep(c(-1,1),each=7))+
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
    geom_text(aes(label = df, size = 5), col = "red")+#, hjust = rep(c(-2,2),each = 7), vjust = rep(c(-1,1),each=7))+
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
ggarrange(gtune[[1]],gtune[[5]],gtune[[2]],gtune[[6]],gtune[[3]],
          gtune[[7]],gtune[[4]],gtune[[8]],nrow = 4, ncol = 2)

##### Final Model & Save Results #####
#Best Hyperparameters
row <- c(6,4,7,7,5,6,8,6)
col <- c(1,1,2,1,1,2,2,2)

set.seed(67489)
#Parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)
registerDoRNG(seed = 9401)

msanet <- list()
msanet_rmse <- list()
msanet_pearson <- list()
test <- list()
metrics <- list()
N <- 100#Bootstrapping
coeff <- list()
for(p in 1:length(y)){
  print(p)
  msanet_rmse[[p]] <- c()
  msanet_pearson[[p]] <- c()
  coeff[[p]] <- data.frame()
  
  features <- read.xlsx(paste0("FS/",names(y)[p],"_Feature_Selection_for.xlsx"), sheetIndex = 1)
  metaboRank$my_var <- y[,p]
  for(i in 1:N){
    print(paste0("N = ", p))
    train_rows <- metaboRank$my_var %>% createDataPartition(p = 0.80, list = F)
    train <- metaboRank[train_rows,c(1,match(features$Feature[features$Counts >= 60], names(metaboRank)))]
    test <- metaboRank[-train_rows,c(1,match(features$Feature[features$Counts >= 60], names(metaboRank)))]
    
    msanet[[p]] <- msaenet(as.matrix(train[,-match("my_var",names(train))]),train[,match("my_var",names(train))],
                           family = "gaussian",init = "enet",tune = "cv",nfolds = 10,
                           rule = criteria[col[p]],alphas = seq(0.1,1,0.1),verbose = T,
                           parallel = T,scale = gamma[row[p]],seed = T,nsteps = steps,tune.nsteps = "max")
    
    #Metrics
    msanet_rmse[[p]] <- c(msanet_rmse[[p]],RMSE(test[,match("my_var",names(test))],predict(msanet[[p]],newx=as.matrix(test[,-match("my_var",names(test))]))))
    msanet_pearson[[p]] <- c(msanet_pearson[[p]],cor.test(test[,match("my_var",names(test))],predict(msanet[[p]],newx=as.matrix(test[,-match("my_var",names(test))])))$est)
    #Coefficients
    cc <- data.frame(names = rownames(msanet[[p]]$beta),coef1 = msanet[[p]]$beta[,1])
    cc <-  cc[cc$coef1 != 0,]
    
    coeff[[p]]  <- merge(coeff[[p]], cc, by = "names", all = T)
    names(coeff[[p]])[length(coeff[[p]])] <- paste0("Iter",i)
  }
  write.xlsx(coeff[[p]], paste0("coef_profiles/Coef_N100_",names(y)[p],".xlsx"))
  write.xlsx(msanet_rmse[[p]], paste0("coef_profiles/Performance_metrics/RMSE_N100_",names(y)[p],".xlsx"))
  write.xlsx(msanet_pearson[[p]], paste0("coef_profiles/Performance_metrics/Pearson_N100_",names(y)[p],".xlsx"))
}
stopCluster(cl)
