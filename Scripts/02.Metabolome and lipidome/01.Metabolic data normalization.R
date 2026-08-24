library(plotly)
library(robustbase)
library(vegan)

#Input metabolomic and background files
data <- read.delim(file.choose(), header = T, row.names = 1)   ##metabolite_matrix.txt
info <- read.delim(file.choose(), header = T, row.names = 1)   ##sample_info.txt

#Permutational multivariate analysis of variance (Before normalization)
permonova<-adonis(data~Group+Gender+Age+Year, data=info, permutations = 9999, method = "euclidean")


#Construct a table for correction coefficients and adjusted results
Adjusted_data = data.frame(matrix(ncol = 362,nrow = nrow(data)))
colnames(Adjusted_data)=colnames(data)
rownames(Adjusted_data)=rownames(data)

Age_Sex_Effects_on_metabolites=data.frame(matrix(ncol = 8,nrow = 362))
colnames(Age_Sex_Effects_on_metabolites)=c("Age_Estimate","Age_SE","Age_T_value","Age_P_value","Sex_Estimate","Sex_SE","Sex_T_value","Sex_P_value")
rownames(Age_Sex_Effects_on_metabolites)=names(data)[1:362]

#Data normalization: linear regression to correct the effects of age and sex
Merge_data <- cbind(data, info, by = "row.names")

for (k in 1:362){
Protein_age_sex_test<-lmrob(Merge_data[,k]~Age+Gender, data=Merge_data, k.max=900000)

Age_Sex_Effects_on_metabolites[k,1:4]=as.vector(summary(Protein_age_sex_test)$coefficients[2,])
Age_Sex_Effects_on_metabolites[k,5:8]=as.vector(summary(Protein_age_sex_test)$coefficients[3,])

Adjusted_data[,k]=Merge_data[,k]-Merge_data[,364]*as.numeric(Age_Sex_Effects_on_metabolites[k,1])-Merge_data[,363]*as.numeric(Age_Sex_Effects_on_metabolites[k,5])
}

#Export correction coefficients and adjusted results
write.csv(Adjusted_data,'Agesex_adjusted_metabolites.csv',row.names = T)
write.csv(Age_Sex_Effects_on_metabolites,'AgeSex_Effects.csv',row.names = T)

#Permutational multivariate analysis of variance (After normalization)

permonova2<-adonis(Adjusted_data~Group+Gender+Age+Year, data=info, permutations = 9999, method = "euclidean")
