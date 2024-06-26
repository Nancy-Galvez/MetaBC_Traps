NODE.MIN <- function (fasta,tree,limite1, limite2, limite3, limite4, limite5, limite6, limite7, limite8, limite9, limite10, limite11, limite12, limiteGMYC,  print.subtrees, print.subtrees.fasta) #un arbol, clase phylo

#Para correr la función usar "defined_thershold.r" dando, en este orden, el fasta y el arbol (nobres identicos) previamente cargados en R, más los 11 limites con el GMYC el último, mas la indicación de YES or NO  para imprimir sibtrees y subfastas  
  
  {
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 1       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  subtrees(tree, wait=TRUE)->subt.t
   
    
     tabla.all.limites <-data.frame()
     
     lista<-data.frame()
    for (i in 1:Nnode(tree))
  {
        if ((max(branching.times(subt.t[[i]])))<limite1)
    {
            for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
        dim(caso)<-c(1,5)
        colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
        lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)

lista_species<-read.table ("lista_species.txt", header=TRUE)
table(lista_species[["spn"]])->nodestab
row.names(nodestab)->nodes_list


#para generar archivo con arbol y fasta para cada subarbol
all.used<-data.frame()
for (h in 1:length(nodes_list))
{
  if (print.subtrees == "YES")
  {
  write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite1, ".nwk") )
  }
  subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
  
  #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
  #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase

subtreecase$tip.label->subtreetips

  #tabla con dos columnas: nombre:seq. sin header
#preparar fasta
as.vector(fasta[,1])->fasta.v

dim(fasta.v)<-c(2,(length(fasta[,1])/2))
t(fasta.v)->fasta.v.t
fasta.v.t->seq_final
seq_final[,1] 

#seq_final
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
#seq_final
GE<-subtreetips  #tabla con lista GE
GE
length(GE)

id_seqs<-data.frame()
for (l in 1:length(GE))
{
  which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
  id_seqs <-c(id_seqs,id)       
}

as.numeric(id_seqs)->id_seqs
id_seqs
#seq_GE<-seq_final[id_seqs,]
#write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 

seq_final_sinGE<-seq_final[id_seqs,]


#trozo para añadir data on a txt
rep(str_c("clade",h,"_limite",limite1),length(seq_final_sinGE[,1]))->used.vector
dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
cbind(seq_final_sinGE,used.vector)->namelist.used

all.used<-rbind(all.used, namelist.used)


#pasarlo a fasta

seq_fasta<-data.frame()
for (g in 1:length(seq_final_sinGE[,1]))
{
  seq_name<-as.vector(seq_final_sinGE[g,1])
  seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
}
dim(seq_fasta)
dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)

#imprimir tabla(fasta format)
if (print.subtrees.fasta == "YES")
{
write.table(seq_fasta, file = str_c("node",h,"_lim",limite1, ".fas"), row.names = FALSE, col.names= FALSE) 
}
#write.table(namelist.used, file=str_c("clade",h,"_limite",limite1,".txt"))
}


#Para generar fasta con las no incluidas en los subtrees
#######################################################################################3
#########################################################################3


id_seqs_unused<-data.frame()
for (k in 1:length(nodes_list))
{
  
  #preparar fasta
  as.vector(fasta[,1])->fasta.v
  
  dim(fasta.v)<-c(2,(length(fasta[,1])/2))
  t(fasta.v)->fasta.v.t
  fasta.v.t->seq_final
  seq_final[,1] 
  
      subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
  #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
  

    id_seqs<-data.frame()
for (y in 1:length(subtreetips))
{
  which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
  id_seqs <-c(id_seqs,id)        
}
id_seqs_unused <-c(id_seqs_unused,id_seqs)
}

as.numeric(id_seqs_unused)->id_seqs_unused
id_seqs
#seq_GE<-seq_final[id_seqs,]
#write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 

seq_final_sinGE<-seq_final[-id_seqs_unused,]

#trozo para añadir data on a txt
rep(str_c("unused_limite", limite1),length(seq_final_sinGE[,1]))->unused.vector
dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
cbind(seq_final_sinGE,unused.vector)->namelist.unused

rbind(all.used, namelist.unused)->all.table.info.clade

#pasarlo a fasta

seq_fasta<-data.frame()
for (g in 1:length(seq_final_sinGE[,1]))
{
  seq_name<-as.vector(seq_final_sinGE[g,1])
  seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
}
dim(seq_fasta)
dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)

#imprimir tabla(fasta format)
if (print.subtrees.fasta == "YES")
{
  write.table(seq_fasta, file = str_c("lim",limite1,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
}
  #write.table(namelist.unused, file=str_c("lim",limite1,"unused", ".txt"))


#datos[order(datosColumna),]


#colnames(all.table.info.clade)<-c("name", "seq", "clade")




write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite1,".txt"))

#para ordenarlos:
read.table (str_c("all.table.info.clade","lim",limite1,".txt"), header=TRUE)->all.table.info.cladeAZ
ordenacion=order(all.table.info.cladeAZ[,1])
all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]

write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite1,".txt"))

tabla.all.limites<-all.table.info.cladeAZ
#}





#################################################################################################################################
#################################################################################################################################
#######################################          LIMITE 2       #################################################################
#################################################################################################################################
#################################################################################################################################

#{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite2)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite2, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite2),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite2, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite2,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite2),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite2,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite2,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite2,".txt"))

  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite2,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite2,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}




#################################################################################################################################
#################################################################################################################################
#######################################          LIMITE 3       #################################################################
#################################################################################################################################
#################################################################################################################################

#{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite3)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite3, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite3),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite3, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite3,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite3),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite3,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite3,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite3,".txt"))

  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite3,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite3,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  #}


#################################################################################################################################
#################################################################################################################################
#######################################          LIMITE 4       #################################################################
#################################################################################################################################
#################################################################################################################################

#{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite4)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite4, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite4),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite4, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite4,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite4),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite4,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite4,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite4,".txt"))

  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite4,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite4,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  #}



#################################################################################################################################
#################################################################################################################################
#######################################          LIMITE 5       #################################################################
#################################################################################################################################
#################################################################################################################################

#{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite5)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite5, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite5),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite5, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite5,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite5),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite5,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite5,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite5,".txt"))

  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite5,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite5,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  #}



  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 6       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite6)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite6, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite6),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite6, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite6,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite6),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite6,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite6,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite6,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite6,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite6,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 7       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite7)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite7, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite7),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite7, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite7,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite7),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite7,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite7,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite7,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite7,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite7,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 8       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite8)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite8, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite8),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite8, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite8,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite8),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite8,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite8,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite8,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite8,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite8,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 9       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite9)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite9, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite9),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite9, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite9,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite9),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite9,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite9,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite9,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite9,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite9,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 10       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite10)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite10, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite10),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite10, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite10,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite10),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite10,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite10,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite10,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite10,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite10,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 11       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite11)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
        dim(caso)<-c(1,5)
        colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
        lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite11, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite11),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite11, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite11,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite11),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite11,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite11,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite11,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite11,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite11,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  
  
  #################################################################################################################################
  #################################################################################################################################
  #######################################          LIMITE 12       #################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  #{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limite12)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
        dim(caso)<-c(1,5)
        colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
        lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limite12, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limite12),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limite12, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limite12,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limite12),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limite12,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limite12,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limite12,".txt"))
  
  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limite12,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limite12,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  
  #}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#################################################################################################################################
#################################################################################################################################
#######################################          LIMITE GMYC       #################################################################
#################################################################################################################################
#################################################################################################################################

#{
  subtrees(tree, wait=TRUE)->subt.t
  
  lista<-data.frame()
  for (i in 1:Nnode(tree))
  {
    if ((max(branching.times(subt.t[[i]])))<limiteGMYC)
    {
      for (x in 1:Ntip(subt.t[[i]]))
      {  
        
        caso<-c("species", i, max(branching.times(subt.t[[i]])),(subt.t[[i]])[[3]][x],Ntip(subt.t[[i]]))
dim(caso)<-c(1,5)
colnames(caso)<-c("sp","spn","depth","taxa", "ntip")
lista<-rbind(lista,caso)
      }
    }
  }
  
  write.table(lista, file="lista.txt", row.names = FALSE, col.names= TRUE)  
  lista<-read.table ("lista.txt", header=TRUE)
  
  tree$tip.label->names
  data.frame()->maxnode
  for (j in 1:length(names))
  {
    subset<-which (lista$taxa == names[j])
    maximum<-max(lista$depth[subset])
    maxtips<-max(lista$ntip[subset])
    asignacion<-which(lista$taxa == names[j] & lista$depth == maximum & lista$ntip == maxtips)
    maxnode<-rbind(maxnode, lista[asignacion,])
  }
  write.table(maxnode, file="lista_species.txt", row.names = FALSE, col.names= TRUE)
  
  lista_species<-read.table ("lista_species.txt", header=TRUE)
  table(lista_species[["spn"]])->nodestab
  row.names(nodestab)->nodes_list
  
  
  #para generar archivo con arbol y fasta para cada subarbol
  all.used<-data.frame()
  for (h in 1:length(nodes_list))
  {
    if (print.subtrees == "YES")
    {
      write.tree ( subt.t[[(as.numeric(nodes_list[h]))]], file = str_c("node",h,"_lim",limiteGMYC, ".nwk") )
    }
    subt.t[[(as.numeric(nodes_list[h]))]]->subtreecase
    
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    
    subtreecase$tip.label->subtreetips
    
    #tabla con dos columnas: nombre:seq. sin header
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    #seq_final
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@-er lista grupos externos y generar fasta sin los GE
    #seq_final
    GE<-subtreetips  #tabla con lista GE
    GE
    length(GE)
    
    id_seqs<-data.frame()
    for (l in 1:length(GE))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(GE[l]))))->id
      id_seqs <-c(id_seqs,id)       
    }
    
    as.numeric(id_seqs)->id_seqs
    id_seqs
    #seq_GE<-seq_final[id_seqs,]
    #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
    
    seq_final_sinGE<-seq_final[id_seqs,]
    
    
    #trozo para añadir data on a txt
    rep(str_c("clade",h,"_limite",limiteGMYC),length(seq_final_sinGE[,1]))->used.vector
    dim(used.vector)<-c(length(seq_final_sinGE[,1]),1)
    cbind(seq_final_sinGE,used.vector)->namelist.used
    
    all.used<-rbind(all.used, namelist.used)
    
    
    #pasarlo a fasta
    
    seq_fasta<-data.frame()
    for (g in 1:length(seq_final_sinGE[,1]))
    {
      seq_name<-as.vector(seq_final_sinGE[g,1])
      seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
    }
    dim(seq_fasta)
    dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
    
    #imprimir tabla(fasta format)
    if (print.subtrees.fasta == "YES")
    {
      write.table(seq_fasta, file = str_c("node",h,"_lim",limiteGMYC, ".fas"), row.names = FALSE, col.names= FALSE) 
    }
    #write.table(namelist.used, file=str_c("clade",h,"_limite",limiteGMYC,".txt"))
  }
  
  
  #Para generar fasta con las no incluidas en los subtrees
  #######################################################################################3
  #########################################################################3
  
  
  id_seqs_unused<-data.frame()
  for (k in 1:length(nodes_list))
  {
    
    #preparar fasta
    as.vector(fasta[,1])->fasta.v
    
    dim(fasta.v)<-c(2,(length(fasta[,1])/2))
    t(fasta.v)->fasta.v.t
    fasta.v.t->seq_final
    seq_final[,1] 
    
    subt.t[[(as.numeric(nodes_list[k]))]]->subtreecase
    #subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
    #subt.t[[(as.numeric(nodes_list[1]))]]->subtreecase
    subtreecase$tip.label->subtreetips
    
    
    id_seqs<-data.frame()
    for (y in 1:length(subtreetips))
    {
      which (seq_final[,1] == as.vector(str_c(">",as.vector(subtreetips[y]))))->id
      id_seqs <-c(id_seqs,id)        
    }
    id_seqs_unused <-c(id_seqs_unused,id_seqs)
  }
  
  as.numeric(id_seqs_unused)->id_seqs_unused
  id_seqs
  #seq_GE<-seq_final[id_seqs,]
  #write.table(seq_GE, file="seq_GE.txt", row.names = FALSE, col.names= FALSE) 
  
  seq_final_sinGE<-seq_final[-id_seqs_unused,]
  
  #trozo para añadir data on a txt
  rep(str_c("unused_limite", limiteGMYC),length(seq_final_sinGE[,1]))->unused.vector
  dim(unused.vector)<-c(length(seq_final_sinGE[,1]),1)
  cbind(seq_final_sinGE,unused.vector)->namelist.unused
  
  rbind(all.used, namelist.unused)->all.table.info.clade
  
  #pasarlo a fasta
  
  seq_fasta<-data.frame()
  for (g in 1:length(seq_final_sinGE[,1]))
  {
    seq_name<-as.vector(seq_final_sinGE[g,1])
    seq_fasta<-c(seq_fasta,seq_name,as.vector(seq_final_sinGE[g,2]))
  }
  dim(seq_fasta)
  dim(seq_fasta)<-c((2*length(seq_final_sinGE[,1])),1)
  
  #imprimir tabla(fasta format)
  if (print.subtrees.fasta == "YES")
  {
    write.table(seq_fasta, file = str_c("lim",limiteGMYC,"unused", ".fas"), row.names = FALSE, col.names= FALSE) 
  }
  #write.table(namelist.unused, file=str_c("lim",limiteGMYC,"unused", ".txt"))
  write.table(all.table.info.clade, file=str_c("all.table.info.clade","lim",limiteGMYC,".txt"))

  
  #para ordenarlos:
  read.table (str_c("all.table.info.clade","lim",limiteGMYC,".txt"), header=TRUE)->all.table.info.cladeAZ
  ordenacion=order(all.table.info.cladeAZ[,1])
  all.table.info.cladeAZ<-all.table.info.cladeAZ[ordenacion,]
  
  write.table(all.table.info.cladeAZ, file=str_c("all.table.info.clade_ordered","lim",limiteGMYC,".txt"))
  
  
  tabla.all.limites<-cbind(tabla.all.limites, all.table.info.cladeAZ[,1], all.table.info.cladeAZ[,3])
  
  write.table( tabla.all.limites, file=str_c(" tabla.all.limites_",limite1,"_",limite2,"_",limite3,"_",limite4,"_",limite5,"_" ,limite6, limite7,"_",limite8,"_",limite9,"_",limite10,"_",limite11,"_",limite12,"_GMYC",limiteGMYC,".txt"))
  }


#subtrees(MIT.s2.coleop, wait=TRUE)->subt.t
#nodes_list[1]
#as.numeric(nodes_list[1])
#subt.t[[12]]
#subt.t[[(nodes_list[1])]]
