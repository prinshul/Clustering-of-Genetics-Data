require(rjags)
require(ggplot2)
require(reshape2)

k<-3  #intialize k to 3
Z<-sample(1:k,172,replace = T) #select values for categorical variable Z randomly. Labels are 1,2,3
frequency=matrix(, nrow = k, ncol = 64) #intialize frequency matrix K x M
data=read.table("8pop_genos.dat")
for(cluster in 1:k)
{
  indicator = as.numeric(rep(cluster,172)==Z)  #find all indexes where Z vector has value equal to cluster
  for(m in 1:64)
  {
    frequency[cluster,m]<-sum(data[m,]*(indicator))/sum(rep(2,172)*(indicator))  # calculate frequency for mth marker for a given cluster
  }
}

gene.data=list(N=ncol(data),M=nrow(data), X=data,k=k)  
gene.model=jags.model("project.model",gene.data,inits=list(Z=Z,Freq=frequency), n.adapt = 1000, n.chains = 2) 
gene.samps=jags.samples(gene.model,variable.names=c("Z","Freq"),1e4)
gene.codasamps=coda.samples(gene.model,variable.names=c("Z","Freq"),1e4)

max.assignment.vect <- vector()
for(i in 1:172)
{
max.assignment=as.numeric(names(which.max(table(gene.samps$Z[i,,1]))))  #find more common assignment of cluster for each individual
max.assignment.vect=c(max.assignment.vect,max.assignment)               #collect assignment for all individuals
}
max.cluster.assignment<-data.frame(individual=rep(1:172),cluster=max.assignment.vect) # data frame having maximum cluster assignemt for each individual 
labels=read.table("8pop_labels.dat")
table(population=labels$V1,cluster=max.assignment.vect)  # cross tabulate the population labels with the maximum cluster assignment 


df=data.frame(matrix(ncol = k+1, nrow = 0)) #data frame to hold average probability for each population in each cluster
colnames(df)<-c("population",1:k)           
pop<-as.vector(unique(labels$V1))      # for each population
for(po in pop)   # for each population
{
 popcols=which(po==labels$V1)  # find individuals belong to a given population
 vect=vector()
 for(index in popcols)         # for each individual
 {
   vect=c(vect,gene.samps$Z[index,,1])  # retrieve thier cluster labels
 }
 rowv<-vector()
 for(cluster in 1:k)    # for each cluster 
 {
   tvect=as.numeric(rep(cluster,1e4*length(popcols))==vect)  # find all indexes with label 'cluster' 
   prob=sum(vect*tvect)/sum(vect)  #sum of all them and divide by total labels for that individual
   rowv<-c(rowv,c(prob))   # keep on combining it for each cluster
 }
 df[nrow(df)+1,] <- c(po,rowv)  #finally add this row <population, prob for k=1, prob for k=2, prob for k=3> to the data frame
}
#boxplot
yvect<-vector()
xvect<-vector()
for(clust in 1:k)
{
  yvect=c(yvect,c(as.numeric(unlist(df[clust+1]))))
  xvect=c(xvect,c(rep(clust,length(unlist(df[clust+1])))))
}
newdf<-data.frame(y=yvect,x=xvect)  #reshape the dataframe
with(newdf, boxplot(y~x,xlab="cluster", ylab="Average Probability",ylim = c(0, 1))) # create boxplots for each cluster and the averge population probabilty 

msd<-vector()
for(iter in 1:1e4)
{
  cpairs=vector()
  for(clust in 1:(k-1))
  {
  for(clust1 in 2:k)
  {
  if(clust1>clust)  # make pairs of clusters
  {
  msd1=gene.samps$Freq[clust,,iter,1]-gene.samps$Freq[clust1,,iter,1]  #select any two pairs and find their difference
  msd1=msd1^2  #square the difference
  msd1=mean(msd1)  #find its mean
  cpairs=c(cpairs,c(msd1))   # store all such means among pairs
  }
  }
  }
  msd=c(msd,c(mean(cpairs))) # finally avergae the mean among the pairs and do this for each iteration
}
plot(msd,t='l',xlab = "iterations",ylab = "average mean squared frequencies")  #trace plot the averge mean squared diff calculated above


logposterior=vector() # find log-posterior upto a constant of proprtionality
N=ncol(data)
M=nrow(data)
alpha=1
beta=1
for(iter in 1:1e4)  # for each iteration
{
  var=0
  for(m in 1:M)  # for each marker
  {
    for(i in 1:N)  # for each individual
    {
      z=gene.samps$Z[i,iter,1]  # get z label value for an individual i in an iteration
      f=gene.samps$Freq[z,m,iter,1]  # get corresponding frequency
      zvect=gene.samps$Z[,iter,1]     # get vector of z at the iter iteration
      vect=as.numeric(rep(z,172)==zvect)  # find all the indexes where zvect has value z
      p=sum(vect*zvect)/sum(zvect)     #probability for categorical distribution for a given z          
      var=var+data[m,i]*log(f)+(2-data[m,i])*log(1-f)+(alpha-1)*log(f)+
        (beta-1)*log(1-f)+log(p)  # sum of log posterior Z and F
      #  constants are thrown away
    }
  }
  logposterior[iter]=var   # do it for all 1e4 samples
}
print(dic.samples(gene.model,1e4))  # find DIC for given k
