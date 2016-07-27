library(DT)
library(shiny)
library(igraph)
library(plotly)
library(rstackdeque)
library(jsonlite)

source("external/graph_utils.R", local = TRUE)
source("external/makenetjson.R", local = TRUE)
source("external/protein_label_dictionary.R",local = TRUE)
conf <- fromJSON("./www/data/config.json")
mp <- getproteinlabeldict()
#head(initial_data)
graph <- build_initial_graph(conf)
communities <- get_communities(graph)
htmlloaded = FALSE
s1 <- rstack()
s2 <-rstack()



function(input, output, session){ 
  global <- reactiveValues()
  global$is_comm_graph = TRUE
  global$viz_stack <- insert_top(s1, list(graph, communities))
  global$name <- insert_top(s2, "")
  
  
  # reset button
  observeEvent(input$reset_button, {
    graph <- build_initial_graph(conf)
    #print(input$select)
    communities <- get_communities(graph,input$select)
    global$viz_stack <- rstack()
    global$viz_stack <- insert_top(global$viz_stack, list(graph, communities))
    global$name <- insert_top(s2, "")
  })
  
  observeEvent(input$variable, {
    #print(input$variable)
  })
  
  #Search button
  observeEvent(input$search_button,{
    searchelm <- strsplit(input$searchentitiy,",")
    #print(searchelm)
    data <- peek_top(global$viz_stack)
    graph <- data[[1]]
    communities <- data[[2]]
    memcomm <- NULL
    if (global$is_comm_graph){
      ii<-1
      for(elm in unlist(searchelm)){
        print(elm)
        memcomm[ii] <-  communities$membership[which(elm== V(graph)$name)]
        ii<-ii+1
      }
      memcommunity<-paste(memcomm,collapse = ",")
    } else {
      memcommunity <- input$searchentitiy
      
    }
   
    observe({
      session$sendCustomMessage(type = "commmemmsg" ,
                                message = list(id=memcommunity))
    })
  })
  
  # table click
  observe({
    row <- input$degree_table_rows_selected
    if (length(row)){
      #print(row)
      session$sendCustomMessage(type = "commmemmsg" ,
                                message = list(id=tail(row, n=1)))
    }
  })
  
  
  # back button
  observeEvent(input$back_button, {
    size <- length(global$viz_stack)
    if (size > 1){
      global$viz_stack <- without_top(global$viz_stack)
      global$name <- without_top(global$name)
    } 
  })
  
  # on-click from sigma.js
  observeEvent(input$comm_id, {
   
	nodeval <- input$comm_id
	matchexp <- paste0("{\"", nodeval,"\":",sep="")
	 con <- file("./www/data/database.json")
	 writesubsr <- NULL
	 open(con)
	 while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
		 print(grepl(matchexp,line,fixed=TRUE))
	 	if((grepl(matchexp,line,fixed=TRUE))==TRUE){
			f <-"./www/data/current_graph.json"
	 		pos = regexpr(matchexp, line,fixed=TRUE)
			writesubsr = substr(line,attr(pos,"match.length") + 1,nchar(line) -1)
			Sys.chmod(f, (file.info(f)$mode | "664"))
			file.remove(f)
			sink(file=f)
			cat(writesubsr)
			sink()
	        
			
	 	}
	 }
	 close(con)
     graph <- build_initial_graph(conf)
     #print(input$select)
     communities <- get_communities(graph,input$select)
     global$viz_stack <- rstack()
     global$viz_stack <- insert_top(global$viz_stack, list(graph, communities))
     global$name <- insert_top(s2, "")
    
  })
  
  # writes out the current viz graph to a json for sigma
  graph_to_write <- reactive({
    data <- peek_top(global$viz_stack)    
    graph <- data[[1]]
    communities <- data[[2]]
    print(global$is_comm_graph)
    # Try and apply community detection if there are a lot of nodes to visualize
    #print(vcount(graph))
    #print(conf$community_threshold)
    if (vcount(graph) >  as.numeric(conf$community_threshold)){
      community_graph <- get_community_graph(graph, communities)
      if (vcount(community_graph) > 1){ 
        global$is_comm_graph <- TRUE
        return(list(community_graph, TRUE))
      }
    } 
    
    # If we have few enough nodes (or would have just 1 (sub)community) visualize as is
    V(graph)$size <- 1
    global$is_comm_graph <- FALSE

    # Remove nodes we aren't we don't want that type of node    
    dellist <- c()
    indx <- 1
    if(input$interactions == "all")
      return(list(graph, FALSE))
    
    for(nd in V(graph)){
      
      atr <- get.vertex.attribute(graph,"type",nd)
      print(atr)
      if(grepl(atr,input$interactions) == FALSE){
        dellist[indx] <- nd
        indx <- indx+1
      }
      
    }
    graph <- delete.vertices(graph,dellist)
    
    return(list(graph, FALSE))
  })
  
  # render with sigma the current graph (in json)
  output$graph_with_sigma <- renderUI({
    data <- graph_to_write()
    makenetjson(data[[1]], "./www/data/current_graph.json", data[[2]],conf) 
    update_stats(data[[1]], data[[2]])
    
    observe({
      session$sendCustomMessage(type = "updategraph",message="xyz")
    })
    
    return(includeHTML("./www/graph.html"))
  })
  
  # update the summary stats
  update_stats <- function(graph, is_comm_graph){
    nodes <- get.data.frame(graph, what="vertices")
    nodes$degree <- degree(graph)
    nodes$pagerank <- page_rank(graph)$vector
    if (is_comm_graph==TRUE){
      colnames(nodes) <- c("Name", "Type", "Comm", "Size", "Degree", "LabelInfo","PageRank")
    } else {
      colnames(nodes) <- c("Name", "Type", "Comm", "Degree", "LabelInfo","PageRank")
    }
    global$nodes <- nodes
  }
  
  # Plot the degree distribution of the current graph
  output$degree_distribution <- renderPlotly({  
    if (!is.null(global$nodes)){
      plot_ly(global$nodes, x = Degree, type="histogram",  color="#FF8800")
    }
  })
  
  # Plot the pagerank distribution of the current graph
  output$pagerank_distribution <- renderPlotly({
    if (!is.null(global$nodes)){
      plot_ly(global$nodes, x = PageRank, type="histogram", color="#FF8800")
    }    
  })
  
  # Generate a table of node degrees
  output$degree_table <- DT::renderDataTable({
    if (!is.null(global$nodes)){
      table <- global$nodes[c("Name", "Degree", "LabelInfo","PageRank")]
    }
  },
  options = list(order = list(list(1, 'desc'))),
  rownames = FALSE,
  selection = "single"
  )
  
  # Generate the current graph name (as a list of community labels)
  output$name <- renderText({
    name <- as.list(rev(global$name))
    name <- paste(name, collapse = "/", sep="/")
    return(paste(c("Current Community", name)))
  })
  
  output$pathway_distribution <- renderPlot({
   # if(global$is_comm_graph == TRUE){
      data <- peek_top(global$viz_stack)
      graph <- data[[1]]
      communities<-data[[2]]
      labellist <- lapply(communities(communities),len)
      rawlabels <- unlist(lapply(labellist,unlist))
      #print(rawlabels)
      labelfreq <- table(rawlabels)
      #lf <- order(labelfreq)[1:10]
      lf <- labelfreq[order(labelfreq,decreasing = T)[1:5]]
      others_cnt <- sum(labelfreq) - sum(lf)
      lf["OTHERS"] <-others_cnt
      pcts <- lapply(lf,function(z){round(100.0*z/sum(lf))})
      pcts <- paste("(",pcts,"%",")",sep="")
      lbls <- paste(pcts,names(lf))
      return(pie(lf,labels = lbls))
      #print(gsize(graph))
      

    #}
    
  })
}
