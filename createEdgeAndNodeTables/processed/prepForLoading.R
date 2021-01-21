node.files <- grep("nodes", dir(), value=TRUE)
edge.files <- grep("edges", dir(), value=TRUE)
tbls.nodes <- lapply(node.files, function(f) read.table(f, sep="\t", as.is=TRUE, header=TRUE))
lapply(tbls.nodes, dim)

tbl.nodes <- unique(do.call(rbind, tbls.nodes))
dim(tbl.nodes)

tbls.edges <- lapply(edge.files, function(f) read.table(f, sep="\t", as.is=TRUE, header=TRUE))
lapply(tbls.edges, dim)

tbl.edges <- unique(do.call(rbind, tbls.edges))
dim(tbl.edges)

filename <- sprintf("nodes-%s.tsv", gsub(" ", ".", date()))
path <- file.path("../../neo4j/import", filename)
write.table(tbl.nodes, file=path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

filename <- sprintf("edges-%s.tsv", gsub(" ", ".", date()))
path <- file.path("../../neo4j/import", filename)
write.table(tbl.edges, file=path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


