
library(igraph)
library(utils)
library(ggplot2)
library(foreach)
library(data.table)
library(dplyr)
library(doParallel)
library(RPostgreSQL)
library(itertools)
library(RColorBrewer)

set.seed(27042012)
nnodes=50
qval_cut=0.1

drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = "yourdb") #db placeholder

d=read.delim("../NEAT_subnetwork_enrichment.phenotypic_factors.tsv",sep="\t",stringsAsFactors=TRUE)
d$a=factor(gsub("\\.(up|down)$","",as.character(d$A),perl=TRUE))
levels(d$a)=unlist(lapply(strsplit(levels(d$a),"\\."),function(x) {if(length(x)==3){return(paste(paste(x[1],x[2],sep="_"),x[3],sep="."))} else{paste(x[1],x[2],sep=".")}}))
D=droplevels(subset(d,conclusion=="Overenrichment" & fdr<0.1))

cl=parallel::makeCluster(nnodes, type="PSOCK")
doParallel::registerDoParallel(cl)

for (combination in levels(D$a)) { 
	x=unlist(strsplit(combination,"\\."))

	analysis=x[1]
	comparison=x[2]
	NEATenriched=as.character(unique(D[D$a==combination,"B"]))
	print(paste(combination,NEATenriched,collapse=","))
	filepre=paste("FDGENEA",analysis,comparison,sep=".")
	
	n=read.delim("../GENIE3.top10_target_ranks.communities_PLM.with_node_centralities.csv",sep=";",stringsAsFactors=FALSE)
	n=as.data.table(n)
	n=setkey(n,"geneid")
	
	a=as.data.table(read.delim("../GENIE3.top10_target_ranks.ranked_regulator_target_interactions.with_communities.correlations_and_directionality.tsv.gz",sep="\t",stringsAsFactors=FALSE))
	a=setkeyv(a,c("regulatoryGene","targetGene"))
	a=a[!regulatoryGene==targetGene]
	
	n$is_regulator=n$geneid %in% unique(a$regulatoryGene)
	n=n[geneid %in% union(a$regulatoryGene,a$targetGene)]
	
	sql="select geneid,name as gene_name, reg_class,superfamily, description from physcogrn.gene_name "
	N=dbGetQuery(con,sql)
	N$gene_name=gsub("Phypa_\\d+,?","",N$gene_name,perl=TRUE,ignore.case =TRUE)
	N$gene_name[N$gene_name ==""]=NA
	n=merge(n,N,by="geneid",all.x=TRUE)
	
	col=read.csv("../subnetwork_colours.csv",stringsAsFactors=FALSE)
	col
	n$color=col[match(n$community,col$network),"colour"]
	
	sql=sprintf("select * from physcogrn.dek1dge_full where analysis='%s' and comparison='%s'",analysis,comparison)
	print(sql)
	P=dbGetQuery(con,sql)
	print(nrow(P))
	
	p=subset(P,geneid %in% n[community %in% NEATenriched][["geneid"]] & qval_lrt<qval_cut)
	print(nrow(p))
	P$direction[P$qval_lrt>=0.1]="unchanged"
	
	n=merge(n,P[,c("geneid","test_stat","b","qval_lrt","direction")],all.x=TRUE)
	n[is.na(n$test_stat),test_stat:=-1]
	n[is.na(n$b),b:=0]
	n[is.na(n$qval_lrt),qval_lrt:=1]
	n[is.na(n$direction),direction:="unchanged"]
	
	a=merge(merge(a,P[,c("geneid","b")],by.x="targetGene",by.y="geneid", all.x=TRUE),P[,c("geneid","b")], by.x="regulatoryGene",by.y="geneid", all.x=TRUE,suffixes=c("_tar","_reg"))
	a[is.na(b_tar),b_tar:=0]
	a[is.na(b_reg),b_reg:=0]
	
	a=a[a$regulatoryGene %in% n$geneid & a$targetGene %in% n$geneid,] 
	A=a[community %in% as.character(as.roman(1:11))]
	ATF=A[regulatoryGene %in% n[reg_class=="TF"][["geneid"]]]
	aTF=a[regulatoryGene %in% n[reg_class=="TF"][["geneid"]]]
	
	g=graph_from_data_frame(a,directed = TRUE,vertices = n)
	G=graph_from_data_frame(A,directed = TRUE,vertices = n)
	gTF=graph_from_data_frame(aTF,directed = TRUE,vertices = n)
	GTF=graph_from_data_frame(ATF,directed = TRUE,vertices = n)
	
	extend=na.omit(unique(unlist(lapply(seq_len(nrow(p)), function(i) {subcomponent(g,p[i,"geneid"],mode="in")$name[1:3]}))))#get genes themselves and the two closests TFs
	fil=unlist(lapply(extend,function(x) {
	    x %in% p[,"geneid"] | length(intersect(p[,"geneid"],subcomponent(g,x,mode="out")$name[1:11]))>1
	}))
	extend=extend[fil]
	
	subG=induced_subgraph(g, extend)
	print(paste(combination, length(V(subG)$name)))
	
	
	extendo=foreach(i=seq_along(extend),.combine="rbind.data.frame",.options.rsr=list(chunkSize=200),.packages=c("igraph"))  %dopar%  {
	    downstream=subcomponent(subG,extend[i],mode="out")$name
	    direct=adjacent_vertices(subG, extend[i], mode = "out")[[1]]$name
	    downstreamab=sum(abs(p[p$geneid %in% downstream, "b"]))
	    downstreamb=sum(p[p$geneid %in% downstream, "b"])
	    directab=sum(abs(p[p$geneid %in% direct, "b"]))
	    directb=sum(p[p$geneid %in% direct, "b"])
	    return(data.frame(
	        geneid=extend[i],
	        associated=extend[i] %in% p[,"geneid"],
	        downstream=length(downstream)-1,
	        downstream_reg=nrow(N[N$geneid %in% downstream & !is.na(N$reg_class), ]),
	        downstream_TF=nrow(N[N$geneid %in% downstream & !is.na(N$reg_class) & N$reg_class=="TF", ]),
	        downstream_cab=downstreamab,
	        downstream_cb=downstreamb,
	        direct=length(direct),
	        direct_reg=nrow(N[N$geneid %in% direct & !is.na(N$reg_class), ]),
	        direct_TF=nrow(N[N$geneid %in% direct & !is.na(N$reg_class) & N$reg_class=="TF", ]),
	        direct_cab=directab,
	        direct_cb=directb
	    ))
	}
	
	extendo$geneid=as.character(extendo$geneid)
	
	extendo=merge(extendo,n,by="geneid")
	
	col2=rep(NA,length(V(subG)))
	col2[V(subG)$is_regulator]="blue"
	col2[V(subG)$is_TF]="black"
	set.seed(27042012)
	lay=layout_with_graphopt(subG)
	pdf(paste(filepre,"full.pdf",sep="."),width=15,height=15)
	plot(subG,layout=lay,vertex.size=2,vertex.frame.color=col2,vertex.label=V(subG)$gene_name,vertex.label.cex=0.5,vertex.color=V(subG)$color,edge.width=0.1,edge.arrow.size=0.1)
	dev.off()
	
	clu=components(subG)
	
	GR=as.data.frame(clu$membership)
	GR=data.frame(geneid=as.character(row.names(GR)),GR)
	row.names(GR)=NULL
	names(GR)[2]="comp"
	
	GR=GR[GR$comp %in% which(table(GR$comp)>1),]
	
	GR$geneid=as.character(GR$geneid)
	GR=merge(GR,extendo,by="geneid")
	
	a=0.5
	b=5
	sizes=a+((log(GR$downstream_cab)-min(log(GR$downstream_cab)))*(b-a)/(max(log(GR$downstream_cab))-min(log(GR$downstream_cab))))
	
	subG2=induced_subgraph(subG, GR[,"geneid"])
	GR=GR[GR$geneid %in% V(subG2)$name,]
	V(subG2)$size<-sizes[match(GR$geneid,V(subG2)$name)]
	
	# Function for 2-step reach
	reach2<-function(x){
	    r=vector(length=vcount(x))
	    for (i in 1:vcount(x)){
	    n=neighborhood(x,2,nodes=i,mode="out")
	    ni=unlist(n)
	    l=length(ni)
	    r[i]=(l)/vcount(x)}
	    names(r)=V(x)$name
	    return(r)
	}
	
	hubscore=hub_score(subG2)$vector
	eigen=evcent(subG2)$vector
	reaching=reach2(subG2)
	outdegree=degree(subG2,mode="out")
	scores=data.frame(reaching,hubscore,eigen,outdegree)
	scores2=data.frame(apply(scores,2,function(x) rank(-x)))
	names(scores2)=paste(names(scores),"rank",sep="_")
	scores=cbind(scores,scores2)
	
	GGR=data.frame(stat=rank(-GR$test_stat),GR)
	
	Q=merge(scores,GGR,by.x=0,by.y="geneid")
	names(Q)[1]="geneid"
	Q=Q[order(-Q$outdegree,-Q$reaching,-Q$hubscore,-Q$eigen,-Q$stat),]
	row.names(Q)=NULL
	
	q=data.frame(unique(Q$community))
	do.call(rbind,lapply(seq_len(nrow(q)),function(i) {
	    pp=Q[Q$community==q[i,1] ,c(19,7,1:6,34,8:13,26,27,28,29)] 
	    if (nrow(pp)>0) { 
	        head(pp,2) }
	}))
	
	set.seed(27042012)
	lay = layout_with_lgl(subG2,root=Q[1,"geneid"],coolexp=30)
	
	COL=V(subG2)$color
	col2=rep(NA,length(V(subG2)))
	col2[V(subG2)$is_TF]="black"
	labs=V(subG2)$gene_name
	
	pdf(paste(filepre,"connected_components.subnetworks.nolab.pdf",sep="."),width=15,height=15)
	plot(subG2,layout=lay,vertex.frame.color=col2,vertex.label=NA,vertex.color=COL,edge.width=0.1,edge.color="black",edge.arrow.size=0)
	dev.off()
	pdf(paste(filepre,"connected_components.subnetworks.lab.pdf",sep="."),width=15,height=15)
	plot(subG2,layout=lay,vertex.frame.color=col2,vertex.label=labs,vertex.label.cex=.7,vertex.color=COL,edge.width=0.1,edge.color="black",edge.arrow.size=0)
	dev.off()
	
	COL=ifelse(V(subG2)$direction=="unchanged", "lightgrey", ifelse(V(subG2)$direction=="up", "red","blue"))
	col2=rep(NA,length(V(subG2)))
	col2[V(subG2)$is_TF]="black"
	labs=V(subG2)$gene_name
	pdf(paste(filepre,"connected_components.direction.nolab.pdf",sep="."),width=15,height=15)
	plot(subG2,layout=lay,vertex.frame.color=col2,vertex.label=NA,vertex.color=COL,edge.width=0.1,edge.color="black",edge.arrow.size=0)
	dev.off()
	pdf(paste(filepre,"connected_components.direction.lab.pdf",sep="."),width=15,height=15)
	plot(subG2,layout=lay,vertex.frame.color=col2,vertex.label=labs,vertex.label.cex=.7,vertex.color=COL,edge.width=0.1,edge.color="black",edge.arrow.size=0)
	dev.off()
	
	write.table(Q,file=paste(filepre,"connected_components.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	write.table(extendo,file=paste(filepre,"all.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	
	options(warn=-1)
	write_graph(subG,file=paste(filepre,"full_graph.xml",sep="."),format="graphml")
	write_graph(subG2,file=paste(filepre,"connected_components_graph.xml",sep="."),format="graphml")
	options(warn=0)
}

stopCluster(cl)
