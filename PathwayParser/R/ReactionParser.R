# ReactionParser.R
#
#' import R6
#' import xml2
#' import EnsDb.Hsapiens.v79
#'
#' @title ReactionParser
#' @description an R6 class which parses xml Reactome reactions into data.frames
#' @name ReactionParser
#'
library(R6)

#' @export
ReactionParser = R6Class("ReactionParser",

    private = list(
        doc = NULL,
        xml.node = NULL,
        parent.pathway = NULL,
        tbl.speciesMap = NULL
        ), # private

        #' @description
        #' Create a new parser
        #' @param doc an xml_document
        #' @param xml.node an XMLInternalNode
        #' @param parent.pathway an object of the Pathway class, constructed from the sbml filename
        #' @return A new `ReactionParser` object.


    public = list(
       initialize = function(doc, xml.node, parent.pathway){
          stopifnot(length(xml.node) == 1)
          private$doc <- doc
          private$xml.node <- xml.node
          private$parent.pathway <- parent.pathway
          private$tbl.speciesMap <- parent.pathway$getMolecularSpeciesMap()
          },


        #' @description
        #' easy access to the entire xml element
        #' @returns an XMLInternalNode
        #' @description
        #' XPath from document root for the reaction
        #' @returns a character string, e.g, "/sbml/model/listOfReactions/reaction[5]"
      getXPath = function(){
         xml_path(private$xml.node)
         },

        #' @description
        #' easy access to the entire xml element
        #' @returns an XMLInternalNode
      getXml = function(){
         private$xml.node
         },

        #' @description
        #' every reaction has a reactome identifier
        #' @returns a character string, the id for this reaction node
      getID = function(){
         xml_attr(private$xml.node, "id")
         },

        #' @description
        #' every reaction has a name
        #' @returns a character string, the name of this reaction node
      getName = function(){
         xml_attr(private$xml.node, "name")
         },

        #' @description
        #' a reaction is specific to a cellular compartmentevery reaction has a name
        #' @returns a character string, the compartment in which  this reaction takes place
      getCompartment = function(){
         xml_attr(private$xml.node, "compartment")
         },

        #' @description
        #' every (most?) reactions are accompanied by explanatory notes
        #' @returns a character string, the short-to-mediums account of the reaction
      getNotes = function(){
         xml_text(xml_find_all(private$xml.node, ".//notes/p"))
         },

        #' @description
        #' one (zero?) or more reactants contribute to each reaction
        #' @returns an integer count
      getReactantCount = function(){
         length(xml_find_all(private$xml.node, ".//listOfReactants/speciesReference"))
         },

        #' @description
        #' get the species identifiers for all reactants
        #' @returns reactant (species) ids
      getReactants = function(){
         xml_text(xml_find_all(private$xml.node, ".//listOfReactants/speciesReference/@species"))
         },

        #' @description
        #' one (zero?) or more products are produced by each reaction
        #' @returns an integer count
      getProductCount = function(){
         length(xml_find_all(private$xml.node, ".//listOfProducts/speciesReference"))
         },

        #' @description
        #' get the species identifiers for all products
        #' @returns product (species) ids
      getProducts = function(){
         xml_text(xml_find_all(private$xml.node, ".//listOfProducts/speciesReference/@species"))
         },

        #' @description
        #' zero or more modifiers are active in each reaction
        #' @returns an integer count
      getModifierCount = function(){
         length(xml_find_all(private$xml.node, ".//listOfModifiers/modifierSpeciesReference"))
         },

        #' @description
        #' get the species identifiers for all modifiers
        #' @returns  character (species) ids
      getModifiers = function(){
         xml_text(xml_find_all(private$xml.node, ".//listOfModifiers/modifierSpeciesReference/@species"))
         },

        #' @description
        #' zero or more molecular complexes participate in this reaction
        #' @returns an integer count
      getComplexCount = function(){
          return(length(self$getComplexes()))
          },

        #' @description
        #' reactants, products, modifiers, complexes
        #' @returns an named list
      getCounts = function(){
          return(list(reactants=self$getReactantCount(),
                      products=self$getProductCount(),
                      modifiers=self$getModifierCount(),
                      complexes=self$getComplexCount()))
          },


        #' @description
        #' identify molecular species which are complexes of molecular species
        #' @returns a named list, complex name and complex participatns
      getComplexes = function(){
         all.species <- unique(c(self$getReactants(), self$getProducts()))
         #tbl.map <- private$parent.pathway$getMolecularSpeciesMap()
         x <- lapply(all.species, function(species) unlist(subset(private$tbl.speciesMap, id==species)$complex.members))
         names(x) <- all.species
         x[as.logical(lapply(x, function(el) !is.null(el)))]
         },


        #' @description
        #' reactome entities (nodes) are of identifiable types; assign them here
        #' @param node.id character
        #' @returns one of "drug", "reaction", "protein", "ligand", or "unrecognized"

      assignNodeType = function(node.id){
          if(node.id == "species_9678687")
              return("drug")
          if(grepl("^species_", node.id))
              return(subset(private$tbl.speciesMap, id==node.id)$type)
          if(grepl("^reaction_", node.id))
              return("reaction")
          if(grepl("^uniprotkb", node.id))
              return("protein")
          if(grepl("^ligandId", node.id))
              return("ligand")
          return("unrecognized")
          }, # assignNodeType

        #' @description
        #' reactome entities (nodes) are of identifiable types; assign them here
        #' @param node.id character
        #' @returns one of "drug", "reaction", "protein", "ligand", or "unrecognized"

      assignNodeName = function(node.id){
          if(node.id %in% c("species_9678687", "ligandId:6031"))
              return("rapamycin")
          if(grepl("^reaction", node.id))
              return(self$getName())
          if(grepl("^uniprotkb", node.id))
              return(select(EnsDb.Hsapiens.v79,
                            key=sub("uniprotkb:", "", node.id),
                            keytype="UNIPROTID",
                            columns=c("SYMBOL"))$SYMBOL)
          if(grepl("^ligandI:", node.id))
              return(node.id)
          if(grepl("^species_", node.id)){
              return(subset(private$tbl.speciesMap, id==node.id)$name)
              }
          return(node.id)
          }, # assignNodeName


        #' @description
        #' cytoscape.js various databases (sql, neo4j, dc) represent data in tables.
        #' create them here
        #' @param excludeUbiquitousSpecies logical, default TRUE: ATP, ADP, water
        #' @param includeComplexMembers logical
        #' @returns a named list, edges and nodes, each a data.frame

      toEdgeAndNodeTables = function(excludeUbiquitousSpecies=TRUE, includeComplexMembers){
          atp <- "species_113592"
          adp <- "species_29370"
          amp <- "species_76577"
          water <- "species_29356"
          inconsequentials <- c(atp, adp, water)
          edge.count <- self$getReactantCount() + self$getProductCount() + self$getModifierCount()
          tbl.in <- data.frame(source=self$getReactants(),
                               target=rep(self$getID(), self$getReactantCount()),
                               interaction=rep("reactantFor", self$getReactantCount()),
                               stringsAsFactors=FALSE)
          tbl.modifiers <- data.frame()
          if(self$getModifierCount() > 0){
              tbl.modifiers <- data.frame(source=self$getModifiers(),
                                         target=rep(self$getID(), self$getModifierCount()),
                                         interaction=rep("modifies", self$getModifierCount()),
                                         stringsAsFactors=FALSE)
             }
          tbl.out <- data.frame(source=rep(self$getID(), self$getProductCount()),
                                target=self$getProducts(),
                                interaction=rep("produces", self$getProductCount()),
                                stringsAsFactors=FALSE)

          tbl.complexes <- data.frame()
          #nodes.complexes <- c()
          if(includeComplexMembers & self$getComplexCount() > 0){
             complex.list <- self$getComplexes()
             tbls.complex <- lapply(names(complex.list), function(species){
                 parent <- species
                 children <- complex.list[[parent]]
                 data.frame(source=children, target=parent, interaction="complexMember",
                            stringsAsFactors=FALSE)
                 })
             tbl.complexes <- do.call(rbind, tbls.complex)
             # nodes.complexes <- c()
             }

          tbl.edges <- rbind(tbl.in, tbl.out, tbl.modifiers) # don't add the complexes
          species <- grep("species_", unique(c(tbl.edges$source, tbl.edges$target)), v=TRUE)
          nodes.all <- with(tbl.edges, unique(c(source, target)))
          if(nrow(tbl.complexes) > 0)
             nodes.all <- c(nodes.all, tbl.complexes$source)
          nodes.species <- intersect(nodes.all, private$tbl.speciesMap$id)
          nodes.other   <- setdiff(nodes.all, private$tbl.speciesMap$id)
          tbl.nodes <- data.frame(id=nodes.all,
                                  type=rep("unassigned", length(nodes.all)), # unlist(lapply(nodes.all, assignNodeType)),
                                  label=unlist(lapply(nodes.all, self$assignNodeName)),
                                  parent=rep("", length(nodes.all)),
                                  stringsAsFactors=FALSE)
          if(nrow(tbl.complexes) > 0){   # add the parents
             complex.members <- match(tbl.complexes$source, tbl.nodes$id)
             tbl.nodes[complex.members, "parent"] <- tbl.complexes$target
             xyz <- 99
             }
              # multiple instances of the same molecule can be distinguished with a .N suffix
              # base R's make.names does this AND converts ":" (e.g., uniprotkb:P62942, ligandId:6031)
              # into period: (e.g., uniprotkb.P62942, ligandId.6031).
              # if used here, those names must be convered back, so that the unmodified nodes
              # mentioned in tbl.edges match up.
          tbl.nodes$id <- make.names(tbl.nodes$id, unique=TRUE)
          tbl.nodes$id <- sub("ChEBI.", "ChEBI:", tbl.nodes$id, fixed=TRUE)
          tbl.nodes$id <- sub("uniprotkb.", "uniprotkb:", tbl.nodes$id, fixed=TRUE)
          tbl.nodes$type <- unlist(lapply(tbl.nodes$id, self$assignNodeType))

          if(excludeUbiquitousSpecies){
             deleters <- match(inconsequentials, tbl.nodes$id)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.nodes <- tbl.nodes[-deleters,]
             deleters <- match(inconsequentials, tbl.edges$source)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.edges <- tbl.edges[-deleters,]
             deleters <- match(inconsequentials, tbl.edges$target)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.edges <- tbl.edges[-deleters,]
             } # if excludeUbiquitousSpecies
          tbl.edges$reaction <- self$getName()
          return(list(edges=tbl.edges, nodes=tbl.nodes))
          }
     ) # public
  ) # class

