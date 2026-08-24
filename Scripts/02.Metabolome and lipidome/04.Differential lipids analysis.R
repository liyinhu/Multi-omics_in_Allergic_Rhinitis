library(ggrepel)
library(ggplot2)

#Import adjusted lipidome
df<-read.csv(file.choose(),header = T)    #adjusted_lipidome.csv

Number_of_lipidites=length(df[1,])-2

Test_summary=data.frame(matrix(ncol=4,nrow=Number_of_lipidites))
colnames(Test_summary)<-c("AR_Estimate","AR_SE","AR_T_value","AR_P_value")
rownames(Test_summary)=names(df[c(3:(Number_of_lipidites+2))])

#Differential analysis between AR and HC using linear regression models
for (i in 1:Number_of_lipidites)
{
tryCatch({
Test<-lm(as.vector(df[,i+2])~Group,data=df)
Test_summary[i,1:4]=as.vector(summary(Test)$coefficients[2,1:4])
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#BH adjustment for multiple testing
Test_summary$AR_FDR=p.adjust(Test_summary$AR_P_value,method='fdr')

#Export statistical results
write.csv(Test_summary, 'Lipidome_summary.csv',row.names = T)

#Volcano plot
k1 <- (Test_summary$AR_P_value < 0.05) &(Test_summary$AR_Estimate > 0)
k2 <- (Test_summary$AR_P_value < 0.05) &(Test_summary$AR_Estimate < 0)
Test_summary$change <- ifelse(k1, "UP",ifelse(k2, "DOWN", "NON"))
Test_summary$label <- ifelse(Test_summary$AR_P_value < 0.05 & (Test_summary$AR_Estimate < -0.575 | Test_summary$AR_Estimate > 0.575), as.character(rownames(Test_summary)), "")

ggplot(data = Test_summary, aes(x = AR_Estimate, y =-log10(AR_P_value), colour =change, fill = change))+ geom_point(alpha =0.6,aes(size=-log10(AR_P_value))+geom_text_repel(aes(x = AR_Estimate, y =-log10(AR_P_value), label =label), size = 3, box.padding = unit(0.6,"lines"),point.padding =unit(0.7,"lines"),  segment.color ="black",show.legend = FALSE)+ geom_hline(yintercept= -log10(0.05),color= 'gray', size = 0.5)+ theme_bw()+ labs(x= "Normalized effect size",y = "-log10(P-value)", title = "Plasma metabolites (AR vs HC)" )+  theme(axis.text = element_text(size = 11), axis.title= element_text(size = 13), plot.title = element_text(hjust = 0.5,size = 15,face ="bold"))+scale_color_manual(values = c('steelblue2','gray30', 'tomato2'))+theme(panel.grid = element_blank(),panel.background = element_rect(color = 'black', fill = 'transparent')) +theme(legend.title = element_blank(),legend.key = element_rect(fill = 'transparent'), legend.background =element_rect(fill = 'transparent'))
