#
#
#
#
rm(list=ls())

library(geotopbricks)
library(dygraphs)


inpts.file <- "geotop.inpts"
##wpath_1D <- '/home/ecor/local/src/geotop_dev/tests/1D' 
wpath_1D <- '../../../tests/1D' 
wpaths <- list.files(wpath_1D,recursive=FALSE,full.name=TRUE)
names(wpaths) <- list.files(wpath_1D,recursive=FALSE,full.name=FALSE)

## Check Simulation Occurence 

inpts.files <- paste(wpaths,inpts.file,sep="/")
igeotopsim <- which(file.exists(inpts.files))

wpaths <- wpaths[igeotopsim]
inpts.files <- inpts.files[igeotopsim]


#wpaths <- wpaths[names(keywords.out)!="Calabria"]
#wpaths <- wpaths[names(keywords.out)!="ColdelaPorte"]
#wpaths <- wpaths[names(keywords.out)!="Jungfraujoch"]




#### Search the output
output_values <- c("output-tabs","output-maps")

keywords <- lapply(X=wpaths,FUN=declared.geotop.inpts.keywords,inpts.file=inpts.file)

timezone <- sapply(X=keywords,FUN=function(x){
			
	        o1 <- str_trim(x$Value[x$Keyword=="MeteoStationStandardTime"])[1]
			o2 <- str_trim(x$Value[x$Keyword=="StandardTimeSimulation"])[1]
			
			print(o1)
			print(o2)
			o <- c(o1,o2)
			o <- o[!is.na(o)][1]
			o <- -as.integer(o)
			if (o>=0) {
				o <- sprintf("Etc/GMT+%d",o)
			} else {
				
				o <- sprintf("Etc/GMT%d",o)
			}
			return(o)
			
		})


keywords.out <- lapply(X=keywords,FUN=function(x,keyword_values){
			cnt <- NULL
			for (kit in keyword_values) {
				
				cnt <- union(cnt,which(str_detect(x$Value,kit)))
				
			}
			
			x[cnt,]
			},keyword_values=output_values)


for (it in names(keywords.out)) {
	
	
	keywords.out[[it]]$sim_name <- it
	keywords.out[[it]]$sim_wpath <- wpaths[it]
	keywords.out[[it]]$formatter <- "%04d"
	keywords.out[[it]]$formatter[keywords.out[[it]]$Keyword %in% c("BasinOutputFile","DischargeFile")] <- ""
	keywords.out[[it]]$tz <- timezone[[it]]
	keywords.out[[it]]$level <- 1
	if (it=="Jungfraujoch") keywords.out[[it]]$level <- "32,33"
}

#keywords.out <- keywords.out[names(keywords.out)!="Calabria"]
#keywords.out <- keywords.out[names(keywords.out)!="ColdelaPorte"]
#keywords.out <- keywords.out[names(keywords.out)!="Jungfraujoch"]
#stop("HERE")
#####
suffixes <- c("-SE27XX","-METEOIO-ON","-METEOIO-OFF")

values.out <- list() 


str(keywords.out)


for (it_s in suffixes) {
	print(it_s)
	values.out[[it_s]] <- lapply(X=keywords.out,FUN=function(x,add_suffix_dir,inpts.file) {

			print(x)
			vars <- x$Keyword
			wpath <- x$sim_wpath[1]
			formatter <- x$formatter
			tz <- x$tz
			
			level <- x$level
			print(level)
			if (is.character(level)) {
				level <- str_split(level,",")
				print(level)
				
			} else {
				
				level <- as.list(level)
			}
			
			print(length(wpath))
			print(length(add_suffix_dir))
			print(length(vars))
			print(length(inpts.file))
			print(length(formatter))
			print(tz)
	#		o <- (mapply(keyword=vars,FUN=get.geotop.inpts.keyword.value,inpts.file=inpts.file,
    #        wpath=wpath,data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",tz=tz,
    #        formatter=formatter,add_suffix_dir=add_suffix_dir,level=level))

			o <- list()
			
			
			for (i in 1:length(vars)) {
				
				for (li in 1:length(level[[i]])) {
				lll <- as.integer(level[[i]][li])
				itn <- sprintf("%s%04d",vars[i],lll)
				o[[itn]] <- get.geotop.inpts.keyword.value(vars[i],inpts.file=inpts.file,
				wpath=wpath,data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",tz=tz[1],
				formatter=formatter[i],add_suffix_dir=add_suffix_dir,level=lll)
				}
			}
			
			
			return(o)
			
			
		},add_suffix_dir=it_s,inpts.file=inpts.file)


}


## Search nemes of the simulations 
sim_names <- names(values.out[[1]])
sim_vars <- NULL
for (it0 in sim_names) {
	
	sim_keywords <- names(values.out[[1]][[it0]])
	for (it1 in sim_keywords) {
		vars <- names(values.out[[1]][[it0]][[it1]])
		sim_vars <- c(sim_vars,paste(it0,it1,vars,sep=";"))
	}
	
}	

sim_vars <<- sim_vars
#print("Testing")
#var <- sim_vars[10]
#
#
#
#data_l <- lapply(X=values.out,FUN=function(x,var){
#			
#			ss <- str_split(var,";")[[1]]
#			sim <- ss[1]
#			kw <- ss[2]
#			col <- ss[3]
#			
#			return(x[[sim]][[kw]][,col])
#			
#			
#			
#		},var=var)
#
#data_time <- index(data_l[[1]])
#names_n <- names(data_l)
#data_l <- lapply(X=data_l,FUN=as.vector)
#data <- as.zoo(as.data.frame(do.call(what="cbind",args=data_l)))
#index(data) <- data_time 
#
#
#
#stop("Testing")





shinyServer(function(input, output) {
#	
#	output$scatterPlot <- renderPlot({
#				x <- rnorm(input$n)
#				y <- rnorm(input$n)
#				plot(x, y)
#			})
	##output$sim_vars <- reactive(c("LOADING",sim_vars))
	output$sim_vars_ui  <-renderUI({selectInput('var', 'Variable',c("LOADING",sim_vars))})
	output$sim_names_ui <-renderUI({selectInput('sim', 'GEOtop Simulation Test',c("LOADING",sim_names))})
	output$sim_kws_ui <-renderUI({selectInput('kws', 'Variable / Keyword',c("LOADING",names(values.out[[1]][[input$sim]])))})
	output$sim_layer_ui <-renderUI({selectInput('layer', 'Layer / Variable',c("LOADING",names(values.out[[1]][[input$sim]][[input$kws]])))})
	
	###uiOutput("dupes")
	
	
	output$dygraph <- 
			renderDygraph({
					sim <- input$sim 
					kws <- input$kws
					layer <- input$layer
					data_l <- lapply(X=values.out,FUN=function(x,sim,kws,layer){x[[sim]][[kws]][,layer]},
								sim=sim,kws=kws,layer=layer)
						
					data_time <- index(data_l[[1]])
					names_n <- names(data_l)
					data_l <- lapply(X=data_l,FUN=as.vector)
					str(data_l)
					data <- as.zoo(as.data.frame(do.call(what="cbind",args=data_l)))
					index(data) <- data_time 
						
					###	residuals <- data-data[,"-SE27XX"]
						
					str(data)
					str(residuals)
					colors <- RColorBrewer::brewer.pal(3, "Set1")
						
						
					dygraph(data, ylab="[unit]") %>% dyRangeSelector() %>%
								dyRoller() %>%
								dyOptions(colors = colors)
						
					})
					
					
			output$dygraph_res <-		
					renderDygraph({
								sim <- input$sim 
								kws <- input$kws
								layer <- input$layer
								data_l <- lapply(X=values.out,FUN=function(x,sim,kws,layer){x[[sim]][[kws]][,layer]},
										sim=sim,kws=kws,layer=layer)
								
								data_time <- index(data_l[[1]])
								names_n <- names(data_l)
								data_l <- lapply(X=data_l,FUN=as.vector)
								str(data_l)
								data <- as.zoo(as.data.frame(do.call(what="cbind",args=data_l)))
								index(data) <- data_time 
								
								###	residuals <- data-data[,"-SE27XX"]
								
								str(data)
								str(residuals)
								
								residuals <- data-data[,"-SE27XX"]
								colors <- RColorBrewer::brewer.pal(3, "Set1")
								dygraph(residuals, ylab="[unit]") %>% dyRangeSelector() %>%
										dyRoller() 		 %>%
										dyOptions(colors = colors)
								
							})
			  output$info <- renderText({
								"bla-bla"
							})
			  
	
})

#
#
#shinyApp(
#		ui = fluidPage(
#				fluidRow(
#					selectInput('var', 'Variable',sim_vars)
#					#	column(4, selectInput('ycol', 'Y Variable', names(dataset),
#					#					selected=names(dataset)[[2]])),
#					#	column(4, numericInput('clusters', 'Cluster count', 3,
#					##					min = 1, max = 9))
#				),
#				fluidRow(
#						plotOutput('kmeans', height = "400px")  
#				)
#		),
#		
#		server = function(input, output, session) {
#			
#			# Combine the selected variables into a new data frame
##			selectedData <- reactive({
##						dataset[, c(input$xcol, input$ycol)]
##					})
##			
##			clusters <- reactive({
##						kmeans(selectedData(), input$clusters)
##					})
##			
##			output$kmeans <- renderPlot(height = 400, {
##						par(mar = c(5.1, 4.1, 0, 1))
##						plot(selectedData(),
##								col = clusters()$cluster,
##								pch = 20, cex = 3)
##						points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
##					})
#		},
#		
#		options = list(height = 500)
#)
#
#stop("CIAO")
#
#

#
#
#
#
#
#
#kmeans_cluster <- function(dataset) { 
#	
#	require(shiny)  
#	
#	shinyApp(
#			ui = fluidPage(
#					fluidRow(style = "padding-bottom: 20px;",
#							column(4, selectInput('xcol', 'X Variable', names(dataset))),
#							column(4, selectInput('ycol', 'Y Variable', names(dataset),
#											selected=names(dataset)[[2]])),
#							column(4, numericInput('clusters', 'Cluster count', 3,
#											min = 1, max = 9))
#					),
#					fluidRow(
#							#plotOutput('kmeans', height = "400px")  
#					)
#			),
#			
#			server = function(input, output, session) {
#				
##				# Combine the selected variables into a new data frame
##				selectedData <- reactive({
##							dataset[, c(input$xcol, input$ycol)]
##						})
##				
##				clusters <- reactive({
##							kmeans(selectedData(), input$clusters)
##						})
##				
##				output$kmeans <- renderPlot(height = 400, {
##							par(mar = c(5.1, 4.1, 0, 1))
##							plot(selectedData(),
##									col = clusters()$cluster,
##									pch = 20, cex = 3)
##							points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
##						})
#			},
#			
#			options = list(height = 500)
#	)
#}
#
#
#
#shinyApp(
#		ui = fluidPage(
#				shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
#				sidebarPanel(
#						textInput("txt", "Text input:", "text here"),
#						selectInput(inputId = "variable", label = "discover variable", choices = unique(c("mela","pera"))),
#						sliderInput("slider", "Slider input:", 1, 100, 30),
#						actionButton("action", "Button"),
#						actionButton("action2", "Button2", class = "btn-primary")
#				),
#				mainPanel(
#						tabsetPanel(
#								tabPanel("Tab 1"),
#								tabPanel("Tab 2")
#						)
#				)
#		),
#		server = function(input, output) {
#			
#			#output$variable <- unique(c(input$variable,input$txt))
#			####ss <<- output$vaqriable
#			####print(output)
#		}
#)
#
#
#
######
##keywords.out <- melt(keywords.out)
##
#
#
#
#
#
#
#
#
