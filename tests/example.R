# adjacency information is in data/germany.adjacency
loc="data/germany.adjacency"

# first number in germany.adjacency is the number of regions
n <- as.numeric(readLines(loc, n = 1))
nnodes <- id <- numeric(n)

# we construct a list of the adjacency information
adj <- list()
# information per line is:
# id of region (id)
# number of neighbours (nnodes)
# ids of the neighbours (adj) 
# attention: ids start with zero!
for (i in 1:n) {
    tmp <- as.numeric(scan(loc, skip = i, nlines = 1, quiet = T, 
                           what = list(rep("", 13)))[[1]])
    id[i] <- tmp[1]
    nnodes[i] <- tmp[2]
    adj[[i]] <- tmp[-c(1:2)]
}

# re-order list per id
adj <- adj[order(id)]
nnodes <- nnodes[order(id)]

# construct a (sparse) matrix from this adjacency structure using sparseMatrix from Matrix package
# define entries which are not zero:
# that is, each entry in the adjacency structure (ids start with zero!)
adj1 <- as.integer(unlist(adj))+1
# plus: diagonal elements are not zero
col     <- c(1:n,  adj1)
# for rows: diagonal elements, id of region
row     <- c(1:n,  rep(1:n,nnodes))
# non-zeor entries: number of neighbours on diagonal, -1 at rest
entries <- c(nnodes, rep(-1,length(adj1)))
# but sparseMatrix only wants upper triagonal matrix
remove <- col<row
entries<-entries[!remove]
col<-col[!remove]
row<-row[!remove]

library(Matrix)
A = sparseMatrix(i=row,j=col,x=entries,symmetric=TRUE) # symmetric=TRUE for symetric matrices

# lets have a look on A
image(A)

# example: random draw
# rmvnorm does draw using cholesky decomposition, see Rue, Held
source("R/rmvnorm.R")
# A has not full rank, we add something to the diagonal:
A1 <- A+0.5*diag(n)
# draw
x <- rmvnorm(1, b=rep(0,n), P=A1)

# germany.plot can be used for plotting the map
library(spam)
germany.plot(x)

# offical ids of Landkreise are given in germany.info$id in spam package
official.ids<-germany.info$id
# id is state id * 1000 + local id
state <- floor(official.ids/1000)

# plot states
colors<-sample(rainbow(16))
germany.plot(state,col=colors,legend=FALSE)

state.names <- read.table("data/states.txt")
legend(1,8800,col=colors,state.names$V2,cex=.5,pch=19)

p <- f(x)
y <- rbinom(544, 1000, p)