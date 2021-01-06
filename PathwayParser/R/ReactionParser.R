# ReactionParser.R
#
#' import R6
#' import xml2
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
        parent.pathway = NULL

        #molecularSpeciesMap = NULL,
        #tbl.pathwayNameMap = NULL,

        # createMolecularSpeciesMap = function(doc){
        #    private$molecularSpeciesMap <- list()
        #    all.species <- xml_find_all(private$doc, "..//listOfSpecies/species")
        #    for(species in all.species){
        #       id <- xml_text(xml_find_all(species, ".//@id"))
        #       name.raw <- xml_text(xml_find_all(species, ".//@name"))
        #       tokens <- strsplit(name.raw, " [", fixed=TRUE)[[1]]
        #       name <- tokens[1]
        #       moleculeType <- "molecule"
        #       compartment <- sub("]", "", tokens[2], fixed=TRUE)
        #       members <- xml_find_all(species, ".//bqbiol:hasPart//rdf:li/@rdf:resource")
        #       if(length(members) == 0)
        #         members <- c()
        #       if(length(members) > 0){
        #         moleculeType <- "complex"
        #         members <- xml_text(members)
        #         members <- sub("http://purl.uniprot.org/uniprot/", "uniprotkb:", members, fixed=TRUE)
        #         members <- sub("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:", "ChEBI:", members, fixed=TRUE)
        #         members <- sub("http://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=",
        #                        "ligandId:", members, fixed=TRUE)
        #         }
        #       new.entry <- list(name=name, moleculeType=moleculeType, compartment=compartment, members=members)
        #       private$molecularSpeciesMap[[id]] <- new.entry
        #       xyz <- 99
        #       } # for species
        #    } # createSpeciesMap
        #readPathwayIdMap = function(){
        #   file <- system.file(package="PathwayParser", "extdata", "ReactomePathways-human.tsv")
        #   tbl.pathwayNameMap <- read.table(file, sep="\t", header=FALSE, as.is=TRUE, nrow=-1, quote="")
        #   colnames(tbl.pathwayNameMap) <- c("id", "name", "organism")
        #   private$tbl.pathwayNameMap <- tbl.pathwayNameMap
        #   }


        ), # private

        #' @description
        #' Create a new parser
        #' @param doc an xml_document
        #' @param xml.node an XMLInternalNode
        #' @return A new `ReactionParser` object.


   public = list(
       initialize = function(doc, xml.node, parent.pathway){
          stopifnot(length(xml.node) == 1)
          private$doc <- doc
          #private$createMolecularSpeciesMap()
          #private$readPathwayIdMap()
          private$xml.node <- xml.node
          private$parent.pathway <- parent.pathway
          },

        # @description
        # associates pathway name with an R-HSA- number
        # @returns the whole data.frame for human pathways
      #getPathwayNameMap = function(){
      #   private$tbl.pathwayNameMap
      #   },

        # @description
        # extracts all molecular species from xml hierarchy into a usable list
        # @returns a named list, indexed by species ids, with reactome data in each element
      #getMolecularSpeciesMap = function(){
      #   private$molecularSpeciesMap
      #   },

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

          tbl.edges <- rbind(tbl.in, tbl.out, tbl.modifiers)

          species <- grep("species_", unique(c(tbl.edges$source, tbl.edges$target)), v=TRUE)
          map <- private$parent.pathway$getMolecularSpeciesMap()

          nodes.all <- with(tbl.edges, unique(c(source, target)))
             # add all the complex members

          nodes.species <- intersect(nodes.all, names(map))
          nodes.other   <- setdiff(nodes.all, names(map))
          assignNodeType <- function(node.id){
              if(node.id == "species_9678687")
                  return("drug")
              if(grepl("^reaction_", node.id))
                  return("reaction")
              if(grepl("^uniprotkb:", node.id))
                  return("protein")
              if(grepl("^ligandId:", node.id))
                  return("ligand")
              if(grepl("^species_", node.id))
                  return(map[[node.id]]$moleculeType)
              return("unrecognized")
          } # assignNodeType

          assignNodeName <- function(node.id){
              if(node.id %in% c("species_9678687", "ligandId:6031"))
                  return("rapamycin")
              if(grepl("^reaction_", node.id))
                  return(self$getName())
              if(grepl("^uniprotkb:", node.id))
                  return(select(EnsDb.Hsapiens.v79,
                                key=sub("uniprotkb:", "", node.id),
                                keytype="UNIPROTID",
                                columns=c("SYMBOL"))$SYMBOL)
              if(grepl("^ligandId:", node.id))
                  return(node.id)
              if(grepl("^species_", node.id)){
                  return(subset(map, id==node.id)$name)
                  }
           return(node.id)
          } # addingNodeName

          tbl.nodes <- data.frame(id=nodes.all,
                                  type=unlist(lapply(nodes.all, assignNodeType)),
                                  label=unlist(lapply(nodes.all, assignNodeName)),
                                  parent=rep("", length(nodes.all)),
                                  stringsAsFactors=FALSE)
          complex.members <- lapply(species, function(sp) map[[sp]]$members)
          names(complex.members) <- species
          tbls.complex <- lapply(names(complex.members), function(parent.node){
              member.proteins <- complex.members[[parent.node]]
              data.frame(id=member.proteins,
                         type=rep("complexMember", length(member.proteins)),
                         label=member.proteins,
                         parent=rep(parent.node, length(member.proteins)),
                         stringsAsFactors=FALSE)
              })
          tbl.complex <- do.call(rbind, tbls.complex)
          tbl.nodes.all <- tbl.nodes
          if(includeComplexMembers)
             tbl.nodes.all <- rbind(tbl.nodes, tbl.complex)

          tbl.nodes.all$id <- make.names(tbl.nodes.all$id, unique=TRUE)

          if(excludeUbiquitousSpecies){
             deleters <- match(inconsequentials, tbl.nodes.all$id)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.nodes.all <- tbl.nodes.all[-deleters,]
             deleters <- match(inconsequentials, tbl.edges$source)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.edges <- tbl.edges[-deleters,]
             deleters <- match(inconsequentials, tbl.edges$target)
             deleters <- deleters[complete.cases(deleters)]
             if(length(deleters) > 0)
                tbl.edges <- tbl.edges[-deleters,]
             } # if excludeUbiquitousSpecies
          return(list(edges=tbl.edges, nodes=tbl.nodes.all))
          }
     ) # public
  ) # class

