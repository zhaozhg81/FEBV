
load("Data/VarComp_L1_prior_1_df_5_G_1000_a_10_b_1.Rdata")
load("Data/VarComp_L1_prior_3_df_5_G_1000_a_10_b_1.Rdata")
load("Data/VarComp_L1_prior_7_df_5_G_1000_a_4_b_1.Rdata")
load("Data/VarComp_L1_prior_6_df_5_G_1000_a_4_b_1.Rdata")

load("Data/VarComp_L1_prior_1_df_5_G_1000_a_6_b_1.Rdata")
load("Data/VarComp_L1_prior_3_df_5_G_1000_a_6_b_1.Rdata")
load("Data/VarComp_L1_prior_7_df_5_G_1000_a_3_b_1.Rdata")
load("Data/VarComp_L1_prior_6_df_5_G_1000_a_3_b_1.Rdata")



TRIM <- 0.05


round(
      log10( cbind( apply(res$tsse1.sSq,2,mean,trim=TRIM), apply(res$tsse1.ljs,2,mean,trim=TRIM) , apply(res$tsse1.opt,2,mean,trim=TRIM) ,
                   apply(res$tsse1.smy,2,mean,trim=TRIM) , apply(res$tsse1.modified.smy,2,mean,trim=TRIM), apply(res$tsse1.vsh,2,mean,trim=TRIM) ,
                   apply(res$tsse1.modified.vsh,2,mean,trim=TRIM) , apply(res$tsse1.rebayes,2,mean,trim=TRIM) , apply(res$tsse1.Feb,2,mean,trim=TRIM)
             )
            ), digits=2
      )
      

round(
      log10( cbind( apply(res$tsse1.sSq,2,mean), apply(res$tsse1.ljs,2,mean) , apply(res$tsse1.opt,2,mean) ,
                   apply(res$tsse1.smy,2,mean) , apply(res$tsse1.modified.smy,2,mean), apply(res$tsse1.vsh,2,mean) ,
                   apply(res$tsse1.modified.vsh,2,mean) , apply(res$tsse1.rebayes,2,mean) , apply(res$tsse1.Feb,2,mean)
             )
            ), digits=2
      )


round(
      log10( cbind( apply(res$tsse1.sSq,2,median), apply(res$tsse1.ljs,2,median) , apply(res$tsse1.opt,2,median) ,
                   apply(res$tsse1.smy,2,median) , apply(res$tsse1.modified.smy,2,median), apply(res$tsse1.vsh,2,median) ,
                   apply(res$tsse1.modified.vsh,2,median) , apply(res$tsse1.rebayes,2,median) , apply(res$tsse1.Feb,2,median)
             )
            ), digits=2
      )
      

load("Data/FinBay_L1_prior_1_G_1000_a_10_b_1.Rdata")
load("Data/FinBay_L1_prior_2_G_1000_a_-2.5_b_1.Rdata")
load("Data/FinBay_L1_prior_3_G_1000_a_10_b_1.Rdata")
load("Data/FinBay_L1_prior_4_G_1000_a_-2.5_b_1.Rdata")

load("Data/FinBay_L1_prior_1_G_1000_a_6_b_1.Rdata")
load("Data/FinBay_L1_prior_2_G_1000_a_-1.5_b_1.Rdata")
load("Data/FinBay_L1_prior_3_G_1000_a_6_b_1.Rdata")
load("Data/FinBay_L1_prior_4_G_1000_a_-1.5_b_1.Rdata")

load("Data/FinBay_L1_prior_1_G_1000_a_10_b_1.Rdata")
load("Data/FinBay_L1_prior_3_G_1000_a_10_b_1.Rdata")
load("Data/FinBay_L1_prior_7_G_1000_a_4_b_1.Rdata")
load("Data/FinBay_L1_prior_6_G_1000_a_4_b_1.Rdata")

load("Data/FinBay_L1_prior_1_G_1000_a_6_b_1.Rdata")
load("Data/FinBay_L1_prior_3_G_1000_a_6_b_1.Rdata")
load("Data/FinBay_L1_prior_7_G_1000_a_3_b_1.Rdata")
load("Data/FinBay_L1_prior_6_G_1000_a_3_b_1.Rdata")


round(
      log10( cbind( mean(res$tsse1.sSq), mean(res$tsse1.ljs), mean(res$tsse1.opt), mean(res$tsse1.smy),
                   mean(res$tsse1.modified.smy), mean(res$tsse1.vsh), mean(res$tsse1.modified.vsh), mean(res$tsse1.rebayes), mean(res$tsse1.Feb)
             )
            ), digits=2
      )
      
