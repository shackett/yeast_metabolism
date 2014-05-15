setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")


library(shiny) # HTML interface
library(gplots) # textPlot
library(ggplot2) # ggplotting
library(gridExtra) # grid of ggplot plots


# find all pathways
# point each pathway to all relevent reactions using a list - order as factors with ALL being first
# generate figures and dump out when appropriate

##### Four files loaded upfront #####
# pathwaySet - number of members in various pathways
# rxToPW - which pathways is a reaction associated with + name
# reactionInfo - reactiion -> rxMechanism + sanitized names
# pathway_plot_list - plots showing information (performance) at a pathway level

#### Reaction-specific plots are loaded on-the-fly when a reaction is chosen #####
# For each reaction (e.g. r_0001), a .Rdata file exists in shinyapp/reaction_data which contains all of the relevent plots

load("shinyapp/shinyData.Rdata", envir=.GlobalEnv)

#test <- F
#if(test == T){
#  load("shinyapp/shinySubData.Rdata") # this is a reduced environment for development purposes
#}else{
#  load("shinyapp/shinyData.Rdata") # this is the full environment
#}



shinyServer(function(input, output) {
  
  if(!file.exists("rOCAplots")){
    dir.create("rOCAplots")
  }
  
  # pathway selected
  pathways_selected <- reactive({input$pathway})
  
  # reactions available
  reactions_available <- reactive({sort(rxToPW$reactionName[rxToPW$pathway == pathwaySet$pathway[pathwaySet$display == pathways_selected()]])})
  output$reactions_available <- renderUI({selectInput("reaction_chosen", "Reaction:", as.list(reactions_available()))})
  
  #viewType = reactive({ifelse(input$reaction_chosen == "PATHWAY INFORMATION", "PW", "RX")})
  
  ### If a pathway is chosen ###
  
  PWplots <- reactive({names(pathway_plot_list[[pathways_selected()]])}) # which pathway is selected
  
  output$pw_check <- renderUI({selectInput("pathway_plots", "Pathway plots to choose", choices = as.list(PWplots()))}) # pass which pathway plots are available to UI
  
  chosenPWplots <- reactive({pathway_plot_list[[pathways_selected()]][names(pathway_plot_list[[pathways_selected()]]) %in% input$pathway_plots]}) # which pathway plots is selected (always just one)
  
  # Update the pathway-level plot that is displayed
  
  output$PW <- renderPlot({
    if(length(chosenPWplots) != 0){
      pwsubList <- chosenPWplots()
      do.call(grid.arrange,  pwsubList) # grid.arrange is used in case future updates result in a variable number of displayed plots
      
    }else{
      return(NULL) 
    }
  })
  
  # If pw_save is triggered, the counter incriments by one and by virtue of changing executes the else statement
  # this saves a ggplot object with a pathway-level summary
  
  observe({
    if(input$pw_save == 0){
      return()
    }else{
      isolate({
        
        pwsubList <- chosenPWplots()
        
        name <- paste0("rOCAplots/", input$pw_filename, ".pdf")
        pdf(file = name, height = 15, width = 15)
        do.call(grid.arrange,  pwsubList)
        dev.off()
      })
    }
  })
  
  
  ### If a reaction is chosen ###
  
  rID = reactive({rxToPW$rID[rxToPW$reactionName == input$reaction_chosen][1]}) # which reaction is active
  
  shiny_flux_data <- reactive({
    if(!is.na(rID())){
      load(paste0("shinyapp/reaction_data/", rID(), "plots.Rdata"))
    }
  })
  shiny_flux_data <- shiny_flux_data()
  
  subtypeChoices <- reactive({reactionInfo$Name[reactionInfo$reaction == rID()]}) # which reaction mechanisms are available
  
  output$subforms_available <- renderUI({selectInput("subform_chosen", "Subtype:", as.list(subtypeChoices()))}) # pass which reaction mechanisms are available to UI
  currentRx <- reactive({reactionInfo$rMech[reactionInfo$reaction == rID()][reactionInfo$Name[reactionInfo$reaction == rID()] == input$subform_chosen]}) # which reaction mechanism is active
  
  # Now that the current reaction is known, load the list containing pre-plotted figures
  #RX_plots_list <- shiny_flux_data()
  #output$indicator <- renderText("shiny_flux_data" %in% ls())
  #output$indicator <- renderText(shiny_flux_data())
  #output$indicator <- renderText(names(RX_plots_list()))
  
  RXplots <- reactive({names(shiny_flux_data()[[currentRx()]]$plotChoices)}) # which plots can be chosen for the reaction form of interest
  
  output$rx_check <- renderUI({checkboxGroupInput("reaction_plots", "Reaction plots to choose", choices = as.list(RXplots()))}) # pass which plots are available to UI
  
  chosenRXplots <- reactive({shiny_flux_data()[[currentRx()]]$plotChoices[names(shiny_flux_data()[[currentRx()]]$plotChoices) %in% input$reaction_plots]}) # a subset of plots to be shown for a reaction form
  
  column_number <- reactive({as.numeric(input$col_num)}) # how many columns should reaction-level figures be displayed in (if only 1 figure is shown it will be 1 column)
  
  # Display information - stoichiometry, pathway, enzymes ... for the reaction of interest
  
  output$RXinfo <- renderPlot({
    textplot(shiny_flux_data()[[currentRx()]]$reactionInfo, cex = 3, valign = "top", halign = "left")
  })
  
  # Update the rxn-level plot(s) shown
  
  output$RX <- renderPlot({
    
    if(length(RXplots) != 0){
      rxsubList <- chosenRXplots()
      rxsubList$ncol = ifelse(length(chosenRXplots()) == 1, 1, column_number())
      
      do.call(grid.arrange,  rxsubList)
      
    }else{
      return(NULL) 
    }
  })
  
  # If rxn_save is triggered, the counter incriments by one and by virtue of changing executes the else statement
  # this saves a ggplot object with a rxn-level summary.  Either a 15x15 plot is generated to look nice as a single plot or 
  # a 25x25 layout of multiple plots is formed for a summary
  
  observe({
    if(input$rxn_save == 0){
      return()
    }else{
      isolate({
        
        rxsubList <- chosenRXplots()
        rxsubList$ncol = ifelse(length(chosenRXplots()) == 1, 1, column_number())
        
        figDimensions <- as.numeric(ifelse(length(chosenRXplots()) == 1, 15, 25))
        
        name <- paste0("rOCAplots/", input$rxn_filename, ".pdf")
        pdf(file = name, height = figDimensions, width = figDimensions)
        
        if(length(chosenRXplots()) != 0){
          do.call(grid.arrange,  rxsubList)
        }else{
          textplot("No plots selected")
        }
        dev.off()
        
        
      })
    }
  })
  
  
  
  
})


