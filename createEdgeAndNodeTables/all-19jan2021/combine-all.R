node.files <- grep("nodes", dir(), value=TRUE)
length(node.files)
edge.files <- grep("edges", dir(), value=TRUE)
length(edge.files)

tbls.nodes <- lapply(node.files, function(f)
    read.table(f, sep="\t", as.is=TRUE, header=TRUE, fill=TRUE, quote=""))
#lapply(tbls.nodes, dim)

atbl.nodes <- unique(do.call(rbind, tbls.nodes))
print(dim(tbl.nodes))  # 35700 4

tbl.nodes$label[32798]
tbl.nodes[32798, "label"] <- gsub('"', '', tbl.nodes[32798, "label"])

tbls.edges <- lapply(edge.files, function(f)
    read.table(f, sep="\t", as.is=TRUE, header=TRUE, fill=TRUE, quote=""))
#lapply(tbls.edges, dim)

tbl.edges <- unique(do.call(rbind, tbls.edges))
print(dim(tbl.edges))  # 50737 4
all(tbl.edges$source %in% tbl.nodes$id)

all.edge.nodes <- unique(c(tbl.edges$source, tbl.edges$target))
length(all.edge.nodes)
missing.nodes <- setdiff(all.edge.nodes, tbl.nodes$id)
length(missing.nodes)

deleter.rows <- c(which(tbl.edges$source %in% missing.nodes),
                  which(tbl.edges$target %in% missing.nodes))
if(length(deleter.rows) > 0)
    tbl.edges <- tbl.edges[-deleter.rows,]

edges.to.fix <- grep("\"Activator", tbl.edges$reaction)
if(length(edges.to.fix) > 0){
    tbl.edges$reaction[edges.to.fix] <-  gsub('"', '', tbl.edges$reaction[edges.to.fix])
    }
dim(tbl.edges)

filename <- sprintf("nodes-%s.tsv", gsub(" ", ".", date()))
path <- file.path("../../neo4j/import", filename)
path
write.table(tbl.nodes, file=path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

filename <- sprintf("edges-%s.tsv", gsub(" ", ".", date()))
path <- file.path("../../neo4j/import", filename)
write.table(tbl.edges, file=path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


