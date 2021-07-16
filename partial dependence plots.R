#############################################################
#################partial dependence plots####################
#############################################################

load("irf20times.RData")
load("Omics_otu.Rdata")
library(iRF)
library(randomForest)

#run iRF
colnames(X_omics) <- gsub("\\W","_",colnames(X_omics))
set.seed(70)
fit_omics <- iRF(x=X_omics[X[[11]]$in_id,], y=Y_omics[X[[11]]$in_id], 
                 xtest=X_omics[-X[[11]]$in_id,], ytest=Y_omics[-X[[11]]$in_id], 
                 n.iter=5,
                 ntree=500,
                 n.core=1,
                 n.bootstrap=30,
                 interactions.return=5,
                 bootstrap.forest=T,
                 rit.param=list(depth=5, ntree=500, nchild=2, 
                                class.id=1, class.cut=NULL))

#Partial dependence plots
#ChoE(20:4)_Streptococcus
pd_omics_crc1 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","OTU7013"), which.class = "CRC")

p1 <- ggplot(pd_omics_crc1, aes(x = pd_omics_crc1[,1], y = pd_omics_crc1[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("Streptococcus")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw() +
  ggtitle("ChoE(20:4)_Streptococcus")
ggsave(p1, filename = "~/Desktop/Partial dependence_Streptococcus.png", 
       width=8, height=6)

#ChoE(20:4)_SM(d18:1/24:1) + SM(d18:2/24:0)
pd_omics_crc2 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","SM_d18_1_24_1____SM_d18_2_24_0_"),
                         which.class = "CRC")

p2 <- ggplot(pd_omics_crc2, aes(x = pd_omics_crc2[,1], y = pd_omics_crc2[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("SM(d18:1/24:1) + SM(d18:2/24:0)")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw()+
  ggtitle("ChoE(20:4)_SM(d18:1/24:1) + SM(d18:2/24:0)") 

ggsave(p2, filename = "~/Desktop/Partial dependence_SM.png", 
       width=8, height=6)

#ChoE(20:4)_SM(42:3)
pd_omics_crc3 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","SM_42_3_"), which.class = "CRC")

p3 <- ggplot(pd_omics_crc3, aes(x = pd_omics_crc3[,1], y = pd_omics_crc3[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("SM(42:3)")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw() +
  ggtitle("ChoE(20:4)_SM(42:3)")
ggsave(p3, filename = "~/Desktop/Partial dependence_SM42.png", 
       width=8, height=6)

#ChoE(20:4)_PE(16:0/18:2)
pd_omics_crc4 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","PE_16_0_18_2_"), which.class = "CRC")

p4 <- ggplot(pd_omics_crc4, aes(x = pd_omics_crc4[,1], y = pd_omics_crc4[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("PE(16:0/18:2)")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw() +
  ggtitle("ChoE(20:4)_PE(16:0/18:2)")
ggsave(p4, filename = "~/Desktop/Partial dependence_PE2.png", 
       width=8, height=6)

#ChoE(20:4)_PE(16:0/18:1)
pd_omics_crc5 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","PE_16_0_18_1_"), which.class = "CRC")

p5 <- ggplot(pd_omics_crc5, aes(x = pd_omics_crc5[,1], y = pd_omics_crc5[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("PE(16:0/18:1)")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw() +
  ggtitle("ChoE(20:4)_PE(16:0/18:1)")
ggsave(p5, filename = "~/Desktop/Partial dependence_PE1.png", 
       width=8, height=6)

#ChoE(20:4)_PC(32:1)
pd_omics_crc6 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","PC_32_1_"), which.class = "CRC")

p6 <- ggplot(pd_omics_crc6, aes(x = pd_omics_crc6[,1], y = pd_omics_crc6[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("PC(32:1)")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw()+
  ggtitle("ChoE(20:4)_PC(32:1)") 
ggsave(p6, filename = "~/Desktop/Partial dependence_PC.png", 
       width=8, height=6)

#ChoE(20:4)_Blautia
pd_omics_crc7 <- partial(fit_omics$rf.list[[5]],train = X_omics[X[[11]]$in_id,], 
                         pred.var = c("ChoE_20_4_","OTU5962"), which.class = "CRC")


p7 <- ggplot(pd_omics_crc7, aes(x = pd_omics_crc7[,1], y = pd_omics_crc7[,2], 
                                z = yhat, fill = yhat))+
  xlab("ChoE(20:4)") + ylab("Blautia")+
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(name = "Centered log odds", palette = "RdYlBu") +
  theme_bw() +
  ggtitle("ChoE(20:4)_Blautia")
ggsave(p7, filename = "~/Desktop/Partial dependence_Blautia.png", 
       width=8, height=6)




